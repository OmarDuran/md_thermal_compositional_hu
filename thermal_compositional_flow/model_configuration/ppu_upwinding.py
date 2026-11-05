"""Phase-potential upwinding (PPU) as a model mixin.

Each phase is upwinded by its own potential flux, on faces and on mortar interfaces,

    Phi_a   = darcy_flux + w(rho_a - rho_hat),
    Psi_a   = interface_darcy_flux + w_perp(rho_a - rho_hat),

and every advective flux is a sum over phases,

    total mass   sum_a <lambda_a>_{Phi_a} Phi_a
    component    sum_a <chi_{xi,a} lambda_a>_{Phi_a} Phi_a
    enthalpy     sum_a <h_a lambda_a>_{Phi_a} Phi_a.

There is no upwinded total mobility and no buoyancy pair term. rho_hat is the
fractionally weighted density in the gravity term of darcy_flux, so it cancels in Phi_a.

Only the total-mass formulation (``fractional_flow=False``) is supported; the
mobility-weighted one puts lambda_T inside the permeability, which PPU does not allow.
Place the mixin before ``FlowModelBase`` in the MRO. PorePy's source is not modified.
"""

from __future__ import annotations

from functools import partial
from typing import Callable, Sequence, cast

import numpy as np
import porepy as pp


class PhasePotentialUpwinding(pp.PorePyModel):
    """Genuine PPU for the total-mass compositional flow model."""

    # ------------------------------------------------------------------ setup checks
    def prepare_simulation(self) -> None:
        if self.params.get("fractional_flow", False):
            raise ValueError("PPU requires fractional_flow=False (total-mass formulation).")
        if not self.params.get("enable_buoyancy_effects", True):
            raise ValueError("PPU requires enable_buoyancy_effects=True (rho_hat in darcy_flux).")
        super().prepare_simulation()
        # Upwind with the initial potentials before the first assembly.
        self.update_buoyancy_driven_fluxes()
        self.rediscretize()

    # ------------------------------------------------------------------ directions
    def ppu_key(self, phase: pp.Phase) -> str:
        return "ppu_" + phase.name

    @pp.ad.cached_method
    def ppu_discretization(self, phase: pp.Phase, subdomains: list[pp.Grid]) -> pp.ad.UpwindAd:
        return pp.ad.UpwindAd(self.ppu_key(phase), subdomains)

    @pp.ad.cached_method
    def ppu_interface_discretization(
        self, phase: pp.Phase, interfaces: list[pp.MortarGrid]
    ) -> pp.ad.UpwindCouplingAd:
        return pp.ad.UpwindCouplingAd(self.ppu_key(phase), interfaces)

    @pp.ad.cached_method
    def phase_potential_flux(self, phase: pp.Phase, subdomains: list[pp.Grid]) -> pp.ad.Operator:
        """Phi_a = darcy_flux + w(rho_a - rho_hat) on the faces.

        On Dirichlet boundary faces the gravity term uses the density of the boundary state,
        as the standalone PPU does, instead of the adjacent cell's.
        """
        rho_hat = self.fractionally_weighted_density(subdomains)
        flux = self.darcy_flux(subdomains) + self.density_driven_flux(
            subdomains, phase.density(subdomains) - rho_hat
        )
        correction = self._boundary_gravity_correction(phase, subdomains)
        if correction is not None:
            flux = flux + correction
        flux.set_name("phase_potential_flux_" + phase.name)
        return flux

    def _boundary_gravity_correction(
        self, phase: pp.Phase, subdomains: list[pp.Grid]
    ) -> pp.ad.Operator | None:
        """c_f rho_a(bc) - w_f(rho_a) on Dirichlet faces, zero elsewhere.

        The gravity flux of face f is w_f(rho) = sum_c a_fc rho_c over the cells of its
        stencil (one cell for TPFA, several for MPFA), and c_f = w_f(1) = sum_c a_fc. The
        correction replaces every rho_c by the boundary state's density, so it holds for
        any flux discretization.
        """
        boundary_grids = self.subdomains_to_boundary_grids(subdomains)
        if not boundary_grids:
            return None
        dirichlet = np.concatenate(
            [self.bc_type_fluid_flux(sd).is_dir.astype(float) for sd in subdomains]
        )
        if not dirichlet.any():
            return None
        projection = pp.ad.BoundaryProjection(self.mdg, subdomains, dim=1)
        rho_bc = projection.boundary_to_subdomain @ phase.density(boundary_grids)
        n_cells = sum(sd.num_cells for sd in subdomains)
        unit = pp.wrap_as_dense_ad_array(1.0, grids=subdomains, size=n_cells)
        c_face = self.density_driven_flux(subdomains, unit)
        w_cell = self.density_driven_flux(subdomains, phase.density(subdomains))
        mask = pp.wrap_as_dense_ad_array(
            dirichlet, grids=subdomains, grid_entity=pp.ad.GridEntity.faces
        )
        return mask * (c_face * rho_bc - w_cell)

    @pp.ad.cached_method
    def interface_phase_potential_flux(
        self, phase: pp.Phase, interfaces: list[pp.MortarGrid]
    ) -> pp.ad.Operator:
        """Psi_a = interface_darcy_flux + w_perp(rho_a - rho_hat) on the mortar cells."""
        subdomains = self.interfaces_to_subdomains(interfaces)
        rho_hat = self.fractionally_weighted_density(subdomains)
        flux = self.interface_darcy_flux(interfaces) + self.interface_density_driven_flux(
            interfaces, phase.density(subdomains) - rho_hat
        )
        flux.set_name("interface_phase_potential_flux_" + phase.name)
        return flux

    # The base model calls these three hooks for the buoyancy directions; with PPU they
    # handle the per-phase potentials instead, so the lag/refresh logic is unchanged.
    def set_buoyancy_discretization_parameters(self) -> None:
        for phase in self.fluid.phases:
            key = self.ppu_key(phase)
            for sd, data in self.mdg.subdomains(return_data=True):
                pp.initialize_data(data, key)
                data[pp.PARAMETERS][key].update(
                    {"darcy_flux": np.zeros(sd.num_faces), "bc": self.bc_type_fluid_flux(sd)}
                )
            for intf, data in self.mdg.interfaces(return_data=True):
                pp.initialize_data(data, key)
                data[pp.PARAMETERS][key].update({"darcy_flux": np.zeros(intf.num_cells)})

    def set_nonlinear_buoyancy_discretization(self) -> None:
        subdomains = self.mdg.subdomains()
        interfaces = self.mdg.interfaces(codim=1)
        for phase in self.fluid.phases:
            self.add_nonlinear_discretization(
                self.ppu_discretization(phase, subdomains).upwind()
            )
            if interfaces:
                self.add_nonlinear_discretization(
                    self.ppu_interface_discretization(phase, interfaces).flux()
                )

    def update_buoyancy_driven_fluxes(self) -> None:
        """Store sign(Phi_a) and sign(Psi_a) as the upwind direction of each phase."""
        subdomains = self.mdg.subdomains()
        interfaces = self.mdg.interfaces(codim=1)
        phases = list(self.fluid.phases)
        ops = [self.phase_potential_flux(p, subdomains) for p in phases]
        if interfaces:
            ops += [self.interface_phase_potential_flux(p, interfaces) for p in phases]
        values = self.equation_system.evaluate(ops)

        sd_offsets = np.cumsum([0] + [sd.num_faces for sd in subdomains])
        intf_offsets = np.cumsum([0] + [intf.num_cells for intf in interfaces])
        for i, phase in enumerate(phases):
            key = self.ppu_key(phase)
            for k, sd in enumerate(subdomains):
                data = self.mdg.subdomain_data(sd)
                data[pp.PARAMETERS][key]["darcy_flux"] = values[i][
                    sd_offsets[k] : sd_offsets[k + 1]
                ]
            for k, intf in enumerate(interfaces):
                data = self.mdg.interface_data(intf)
                data[pp.PARAMETERS][key]["darcy_flux"] = values[len(phases) + i][
                    intf_offsets[k] : intf_offsets[k + 1]
                ]

    # ------------------------------------------------------------------ phase sums
    def _ppu_face_flux(
        self,
        subdomains: list[pp.Grid],
        weights: dict[pp.Phase, Callable[[Sequence], pp.ad.Operator]],
        neumann_bc: pp.ad.Operator,
        bc_type: Callable[[pp.Grid], pp.BoundaryCondition],
        neumann_operator: Callable[[Sequence[pp.BoundaryGrid]], pp.ad.Operator],
        interface_flux: Callable[[list[pp.MortarGrid]], pp.ad.Operator],
        name: str,
    ) -> pp.ad.Operator:
        """sum_a <w_a>_{Phi_a} Phi_a with boundary and interface terms.

        Dirichlet faces take the upstream boundary value of w_a along Phi_a; Neumann data
        and the interface flux are total quantities and enter once.
        """
        terms = []
        for phase, weight in weights.items():
            phi = self.phase_potential_flux(phase, subdomains)
            discr = self.ppu_discretization(phase, subdomains)
            bc_phase = self._combine_boundary_operators(  # type: ignore[call-arg]
                subdomains=subdomains,
                dirichlet_operator=weight,
                neumann_operator=neumann_operator,
                robin_operator=None,
                bc_type=bc_type,
                name=f"bc_values_{name}_{phase.name}",
            )
            terms.append(
                phi * (discr.upwind() @ weight(subdomains))
                + discr.bound_transport_dir() @ (phi * bc_phase)
            )
        neu = self.mobility_discretization(subdomains).bound_transport_neu()
        flux = pp.ad.sum_operator_list(terms) + neu @ neumann_bc
        interfaces = self.subdomains_to_interfaces(subdomains, [1])
        if interfaces:
            projection = pp.ad.MortarProjections(self.mdg, subdomains, interfaces, dim=1)
            flux = flux + neu @ projection.mortar_to_primary_int() @ interface_flux(interfaces)
        flux.set_name(name)
        return flux

    def _ppu_interface_flux(
        self,
        interfaces: list[pp.MortarGrid],
        weights: dict[pp.Phase, Callable[[Sequence], pp.ad.Operator]],
        name: str,
    ) -> pp.ad.Operator:
        """sum_a Psi_a <w_a>_{Psi_a} on the mortar cells."""
        subdomains = self.interfaces_to_subdomains(interfaces)
        projection = pp.ad.MortarProjections(self.mdg, subdomains, interfaces, dim=1)
        trace = pp.ad.Trace(subdomains)
        terms = []
        for phase, weight in weights.items():
            discr = self.ppu_interface_discretization(phase, interfaces)
            w = weight(subdomains)
            upstream = (
                discr.upwind_primary() @ projection.primary_to_mortar_avg() @ trace.trace @ w
                + discr.upwind_secondary() @ projection.secondary_to_mortar_avg() @ w
            )
            terms.append(self.interface_phase_potential_flux(phase, interfaces) * upstream)
        flux = pp.ad.sum_operator_list(terms)
        flux.set_name(name)
        return flux

    # Per-phase advected weights: lambda_a, chi_{xi,a} lambda_a, h_a lambda_a.
    def _mass_weights(self) -> dict:
        return {p: partial(self._phase_mass_mobility, p) for p in self.fluid.phases}

    def _component_weights(self, component: pp.Component) -> dict:
        def weight(phase, domains):
            return phase.partial_fraction_of[component](domains) * self._phase_mass_mobility(
                phase, domains
            )

        return {
            p: partial(weight, p) for p in self.fluid.phases if component in p  # type: ignore[operator]
        }

    def _enthalpy_weights(self) -> dict:
        def weight(phase, domains):
            return phase.specific_enthalpy(domains) * self._phase_mass_mobility(phase, domains)

        return {p: partial(weight, p) for p in self.fluid.phases}

    # ------------------------------------------------------------------ fluxes
    @pp.ad.cached_method
    def fluid_flux(self, domains: pp.SubdomainsOrBoundaries) -> pp.ad.Operator:
        if len(domains) == 0 or all(isinstance(d, pp.BoundaryGrid) for d in domains):
            return super().fluid_flux(domains)
        subdomains = cast(list[pp.Grid], list(domains))
        return self._ppu_face_flux(
            subdomains,
            self._mass_weights(),
            self.boundary_fluid_flux(subdomains),
            self.bc_type_fluid_flux,
            self.fluid_flux,
            self.interface_fluid_flux,
            "ppu_fluid_flux",
        )

    @pp.ad.cached_method
    def interface_fluid_flux(self, interfaces: list[pp.MortarGrid]) -> pp.ad.Operator:
        return self._ppu_interface_flux(interfaces, self._mass_weights(), "ppu_interface_fluid_flux")

    def component_flux(
        self, component: pp.Component, domains: pp.SubdomainsOrBoundaries
    ) -> pp.ad.Operator:
        if len(domains) == 0 or all(isinstance(d, pp.BoundaryGrid) for d in domains):
            return super().component_flux(component, domains)
        subdomains = cast(list[pp.Grid], list(domains))
        return self._ppu_face_flux(
            subdomains,
            self._component_weights(component),
            self.boundary_component_flux(component, subdomains),
            self.bc_type_fluid_flux,
            partial(self.component_flux, component),
            partial(self.interface_component_flux, component),
            f"ppu_component_flux_{component.name}",
        )

    def interface_component_flux(
        self, component: pp.Component, interfaces: list[pp.MortarGrid]
    ) -> pp.ad.Operator:
        return self._ppu_interface_flux(
            interfaces,
            self._component_weights(component),
            f"ppu_interface_component_flux_{component.name}",
        )

    def enthalpy_flux(self, domains: pp.SubdomainsOrBoundaries) -> pp.ad.Operator:
        if len(domains) == 0 or all(isinstance(d, pp.BoundaryGrid) for d in domains):
            return super().enthalpy_flux(domains)
        subdomains = cast(list[pp.Grid], list(domains))
        neumann_bc = self._combine_boundary_operators(  # type: ignore[call-arg]
            subdomains=subdomains,
            dirichlet_operator=self.advection_weight_energy_balance,
            neumann_operator=self.enthalpy_flux,
            robin_operator=None,
            bc_type=self.bc_type_enthalpy_flux,
            name="bc_values_enthalpy",
        )
        return self._ppu_face_flux(
            subdomains,
            self._enthalpy_weights(),
            neumann_bc,
            self.bc_type_enthalpy_flux,
            self.enthalpy_flux,
            self.interface_enthalpy_flux,
            "ppu_enthalpy_flux",
        )

    def interface_enthalpy_flux_equation(
        self, interfaces: list[pp.MortarGrid]
    ) -> pp.ad.Operator:
        flux = self._ppu_interface_flux(
            interfaces, self._enthalpy_weights(), "ppu_interface_enthalpy_flux"
        )
        eq = self.interface_enthalpy_flux(interfaces) - flux
        eq.set_name("interface_enthalpy_flux_equation")
        return eq

    # ------------------------------------------------------------------ no pair terms
    def _zero_faces(self, domains: Sequence[pp.Grid]) -> pp.ad.Operator:
        return pp.wrap_as_dense_ad_array(
            0.0,
            grids=list(domains),
            size=sum(sd.num_faces for sd in domains),
            grid_entity=pp.ad.GridEntity.faces,
        )

    def _zero_cells(self, domains: Sequence[pp.Grid]) -> pp.ad.Operator:
        return pp.wrap_as_dense_ad_array(
            0.0, grids=list(domains), size=sum(sd.num_cells for sd in domains)
        )

    def component_buoyancy(self, component, domains) -> pp.ad.Operator:
        return self._zero_faces(domains)

    def enthalpy_buoyancy(self, domains) -> pp.ad.Operator:
        return self._zero_faces(domains)

    def component_buoyancy_jump(self, component, domains) -> pp.ad.Operator:
        return self._zero_cells(domains)

    def enthalpy_buoyancy_jump(self, domains) -> pp.ad.Operator:
        return self._zero_cells(domains)
