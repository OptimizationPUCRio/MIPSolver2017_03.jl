using JuMP, Gurobi
include("functions.jl")
############################### ###############################
############################### ###############################
############################### ###############################
############################### ###############################
xlb = [1,0,1]
xup = [3,1,+Inf]
solver= GurobiSolver(OutputFlag=0)
m = Model(solver=solver)
@variable(m, xlb[i] <= x[i=1:3] <= xup[i], Int)
@constraint(m, 5*x[1] - 2*x[2] + 8*x[3]<= 15)
@constraint(m, 8*x[1] + 3*x[2] - x[3] >= 9)
@constraint(m,x[1] + x[2] + x[3] <= 6)
@objective(m,Max, 2x[1] + x[2] - x[3])
############################### ###############################
############################### ###############################

SolverBrito.SolveMIP(m)

#solve(m)

println("Z = ", getobjectivevalue(m))
println("Vars = ", getvalue(x))
println("Status = ", m.ext[:status])
println("Tempo = ", m.ext[:time])
