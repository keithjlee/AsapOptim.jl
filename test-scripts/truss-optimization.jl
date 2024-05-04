using AsapOptim, Asap, AsapToolkit
using Zygote
using LinearSolve, LinearAlgebra

begin
    w_section = W("W460X158")
    @show w_section.name

    section = Section(
        Steel_kNm,
        w_section.A / 1e6,
        w_section.Ix / 1e12,
        w_section.Iy / 1e12,
        w_section.J / 1e12
    )
end;

# generate
begin
    nx = 25
    dx = 1.5
    ny = 25
    dy = 1.75
    dz = 2.

    # loads
    load = [0., 0., -20]

    spaceframe = SpaceFrame(nx, dx, ny, dy, dz, section; load = load)
    model = spaceframe.model
end

vars = [SpatialVariable(node, 0., -.25dz, dz, :Z) for node in model.nodes[:top]]

params = TrussOptParams(model, vars)