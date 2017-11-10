using JuMP, Gurobi, Base.Test, MathProgBase
include("functions.jl")
############################### ###############################
############################### ###############################
############################### ###############################
############################### ###############################
xlb = [0,0]
xup = [Inf,Inf]
solver= GurobiSolver(OutputFlag=0)
m = Model(solver=solver)
@variable(m, x[i=1:2] >= 0)
@constraint(m, 2x[1] + x[2] <= 4)
@constraint(m, x[1] + 2x[2] <= 4)

@objective(m,Max, 4x[1] + 3x[2])
############################### ###############################
############################### ###############################

SolverBrito.SolveMIP(m)

#solve(m)

println("Z = ", getobjectivevalue(m))
#println("Vars = ", cat((2,2),transpose(getvalue(x)),transpose(getvalue(y))))
println("Status = ", m.ext[:status])
println("Tempo = ", m.ext[:time])
