function Rframe(Cx::Float64, Cy::Float64, Cz::Float64, Ψ::Float64; tol = 1e-4)

    cΨ = cos(Ψ)
    sΨ = sin(Ψ)


    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. Cy 0.;
            -Cy*cΨ 0 sΨ;
            Cy*sΨ 0 cΨ]
    else # all other
        b1 = (-Cx * Cy * cΨ - Cz * sΨ) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cΨ
        b3 = (-Cy * Cz * cΨ + Cx * sΨ) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sΨ - Cz * cΨ) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sΨ
        c3 = (Cy * Cz * sΨ + Cx * cΨ) / sqrt(Cx^2 + Cz^2)

        Λ = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end

function Rframe(Cxyz::SubArray, Ψ::Float64; tol = 1e-4)

    Cx, Cy, Cz = Cxyz

    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. Cy 0.;
            -Cy*cΨ 0 sΨ;
            Cy*sΨ 0 cΨ]
    else # all other
        b1 = (-Cx * Cy * cΨ - Cz * sΨ) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cΨ
        b3 = (-Cy * Cz * cΨ + Cx * sΨ) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sΨ - Cz * cΨ) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sΨ
        c3 = (Cy * Cz * sΨ + Cx * cΨ) / sqrt(Cx^2 + Cz^2)

        Λ = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end

Rframe(XYZn::Matrix{Float64}, Ψ::Float64; tol = 1e-4) = [Rtruss(xyzn, psi; tol = tol) for (xyzn, psi) in zip(eachrow(XYZn), Ψ)]


function kframe(E::Float64, A::Float64, L::Float64, G::Float64, Izz::Float64, Iyy::Float64, J::Float64)
    E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Izz 0 0 0 6L*Izz 0 -12Izz 0 0 0 6L*Izz;
        0 0 12Iyy 0 -6L*Iyy 0 0 0 -12Iyy 0 -6L*Iyy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 4L^2*Iyy 0 0 0 6L*Iyy 0 2L^2*Iyy 0;
        0 6L*Izz 0 0 0 4L^2*Izz 0 -6L*Izz 0 0 0 2L^2*Izz;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Izz 0 0 0 -6L*Izz 0 12Izz 0 0 0 -6L*Izz;
        0 0 -12Iyy 0 6L*Iyy 0 0 0 12Iyy 0 6L*Iyy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 2L^2*Iyy 0 0 0 6L*Iyy 0 4L^2*Iyy 0;
        0 6L*Izz 0 0 0 2L^2*Izz 0 -6L*Izz 0 0 0 4L^2*Izz
    ]
end




# function kframe_01(E::Float64, A::Float64, L::Float64, G::Float64, Izz::Float64, Iyy::Float64, J::Float64)

#     E / L^3 .* [
#         A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
#         0 3Izz 0 0 0 0 0 -3Izz 0 0 0 3L*Izz;
#         0 0 3Iyy 0 0 0 0 0 -3Iyy 0 -3L*Iyy 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
#         0 -3Izz 0 0 0 0 0 3Izz 0 0 0 -3L*Izz;
#         0 0 -3Iyy 0 0 0 0 0 3Iyy 0 3L*Iyy 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 -3L*Iyy 0 0 0 0 0 3L*Iyy 0 3L^2*Iyy 0;
#         0 3L*Izz 0 0 0 0 0 -3L*Izz 0 0 0 3L^2*Izz    
#     ]
# end

# function kframe_10(E::Float64, A::Float64, L::Float64, G::Float64, Izz::Float64, Iyy::Float64, J::Float64)

#     E / L^3 .* [
#         A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
#         0 3Izz 0 0 0 3L*Izz 0 -3Izz 0 0 0 0;
#         0 0 3Iyy 0 -3L*Iyy 0 0 0 -3Iyy 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 -3L*Iyy 0 3L^2*Iyy 0 0 0 3L*Iyy 0 0 0;
#         0 3L*Izz 0 0 0 3L^2*Izz 0 -3L*Izz 0 0 0 0 ;
#         -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
#         0 -3Izz 0 0 0 -3L*Izz 0 3Izz 0 0 0 0;
#         0 0 -3Iyy 0 3L*Iyy 0 0 0 3Iyy 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0    
#     ]
# end

# function kframe_00(E::Float64, A::Float64, L::Float64, G::Float64, Izz::Float64, Iyy::Float64, J::Float64)
#     E * A / L .* [
#         1. 0 0 0 0 0 -1. 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         -1. 0 0 0 0 0 1. 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0;
#         0 0 0 0 0 0 0 0 0 0 0 0
#         ]
# end