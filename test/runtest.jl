using Test
using DynaMeth

@testset "Reactor construction with different kinetics" begin
    for kin in (K() for K in subtypes(AbstractKineticSetup))
        fpo = SetupFpoPftr(
            n=50, ordn=2, m=2, ordm=2, numP=40, N2=0.15,
            kinetic=kin,
            optimize_temperature=false,
            isothermal=true,
        )

        reac = Reactor(fpo)

        @test isa(reac, Reactor)
    end
end;
