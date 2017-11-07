using JuMP, Gurobi, Base.Test, MathProgBase
include("functions.jl")
############################### ###############################
############################### ###############################
############################### ###############################
############################### ###############################
xlb = [0,0,0]
xup = [3,1,+Inf]
solver= GurobiSolver(OutputFlag=0)
m = Model(solver=solver)
@variable(m, xlb[i] <= x[i=1:1] <= xup[i])
@variable(m, xlb[i+1] <= y[i=1:2] <= xup[i+1],Bin)
@constraint(m, 5*x[1] - 2*y[1] + 8*y[2]<= 15)
@constraint(m, 8*x[1] + 3*y[1] - y[2] >= 9)
@constraint(m,x[1] + y[1] + y[2] <= 6)
@objective(m,Max, 2x[1] + y[1] - y[2])
############################### ###############################
############################### ###############################

SolverBrito.SolveMIP(m)

#solve(m)

println("Z = ", getobjectivevalue(m))
println("Vars = ", cat((2,2),transpose(getvalue(x)),transpose(getvalue(y))))
println("Status = ", m.ext[:status])
println("Tempo = ", m.ext[:time])
