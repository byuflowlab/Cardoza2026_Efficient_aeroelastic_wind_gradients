using UnsteadyOpt, Test, StaticArrays, SparseArrays

uo = UnsteadyOpt

# @testset "utils.jl" begin

#     distro = ()
#     assembly = ()
#     aero = ()
#     structural = ()
#     composites = ()
#     idxs = ()
#     scaling = ()

#     obj = uo.ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling)

#     nx = 10
#     ng = 5
#     deriv = uo.Deriv(obj, (nx, ng))
#     # @code_warntype uo.Deriv(obj, (nx, ng))
#     @test isa(@inferred(uo.Deriv(obj, (nx, ng))), uo.Deriv)

#     function (obj::uo.ObjectiveFunction)(g, x)
#         for i in eachindex(g)
#             g[i] = 2*i  # Dummy constraint function
#         end

#         return sum(x) # Dummy objective function
#     end

#     x = rand(nx)
#     g = zeros(ng)
#     f = obj(g, x)
#     # @code_warntype obj(g, x)
#     @test isa(@inferred(obj(g, x)), Float64)

#     df = zeros(nx)
#     dg = zeros(ng, nx)
#     f = deriv(g, df, dg, x)
#     # @code_warntype deriv(g, df, dg, x)
#     @test isa(@inferred(deriv(g, df, dg, x)), Float64)


#     testlist = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
#     list = uo.get_list_property(testlist)
#     # @code_warntype uo.get_list_property(testlist)
#     @test isa(@inferred(uo.get_list_property(testlist)), Vector{Int})
#     @test list == [1, 2, 3, 4, 5, 6, 7, 8, 9]

#     inlist = collect(1:36)
#     properties = uo.create_list_property(inlist) #Todo. Is this the correct order? -> The matrices are symmetric, so if it's column-major, it should be fine or row-major, it doesn't matter. 
#     # @code_warntype uo.create_list_property(inlist)
#     @test isa(@inferred(uo.create_list_property(inlist)), Vector{SMatrix{6, 6, Int, 36}})



#     distro = ()
#     assembly = ()
#     aero = ()
#     structural = ()
#     composites = ()
#     idxs = ()
#     scaling = ()

#     obj = uo.ObjectiveFunction(distro, assembly, aero, structural, composites, idxs, scaling)
#     con = uo.ConstraintFunction(distro, assembly, aero, structural, composites, idxs, scaling)

#     function (obj::uo.ObjectiveFunction)(x)
#         println("called obj")
#         f = x[1]^2 - x[2]

#         return f
#     end

#     function (con::uo.ConstraintFunction)(g, x)
#         # g[1] = x[2] - 2*x[1]
#         println("called con")
#         g[1] = - 2*x[1] #Removing x[2] to test sparsity
#         g[2] = -x[2]
#         g[3] = x[1]^2

#     end

#     nx = 2 
#     ng = 3

#     println("initializing sparse deriv")
#     deriv = uo.SparseDeriv(obj, con, (nx, ng))
    
#     println("Calling objective...")
#     x = [1.0, 2.0]
#     f1 = obj(x)

#     println("Calling constraint...")
#     g = zeros(ng)
#     con(g, x)

#     println("calling deriv")
#     df = zeros(nx)
#     dg = zeros(ng, nx)
#     f3 = deriv(g, df, dg, x)


# # end

@testset "average_every_n tests" begin
    # Test basic functionality
    @test uo.average_every_n([1., 2, 3, 4, 5, 6], 2) == [1.5, 3.5, 5.5]
    @test uo.average_every_n([1., 2, 3, 4, 5, 6], 3) == [2.0, 5.0]
    
    # Test when length is not divisible by n
    @test uo.average_every_n([1., 2, 3, 4, 5], 2) == [1.5, 3.5, 5.0]
    @test uo.average_every_n([1., 2, 3, 4, 5, 6, 7], 3) == [2.0, 5.0, 7.0]
    
    # Test edge cases
    @test uo.average_every_n([5.], 1) == [5.0]  # Single element
    @test uo.average_every_n([1., 2], 1) == [1.0, 2.0]  # n=1 (no averaging)
    @test uo.average_every_n([1., 2., 3., 4.], 4) == [2.5]  # n equals length
    @test uo.average_every_n([1., 2., 3.], 5) == [2.0]  # n > length

    # Test with floating point numbers
    @test uo.average_every_n([1.1, 2.2, 3.3, 4.4], 2) ≈ [1.65, 3.85]
    
    # Test type preservation
    x_float = [1.0, 2.0, 3.0, 4.0]
    result_float = uo.average_every_n(x_float, 2)
    @test eltype(result_float) == Float64
    
    # Test with different vector types
    @test uo.average_every_n(1.:6, 2) == [1.5, 3.5, 5.5]
    @test uo.average_every_n(Vector{Float32}([1, 2, 3, 4]), 2) ≈ [1.5, 3.5]
end



@testset "sparsity" begin
    A = sparse([1 0 3; 0 5 0; 7 0 9])
    B = sparse([1 2 0; 0 5 0; 0 0 9])

    only_in_A, only_in_B = uo.compare_sparsity(A, B)

    @test only_in_A == [(3, 1), (1, 3)]
    @test only_in_B == [(1, 2)]

end