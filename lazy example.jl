using CPLEX
using JuMP

model = JuMP.Model(CPLEX.Optimizer)
@variable(model, x[1:2] <= 10, Int)
@variable(model, y <= 10, Int)
@objective(model, Max, sum(x)+y)


function my_callback_function(cb_data)

    status = callback_node_status(cb_data, model)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        # `callback_value(cb_data, x)` is not integer (to some tolerance).
        # If, for example, your lazy constraint generator requires an
        # integer-feasible primal solution, you can add a `return` here.
        #return
    elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # `callback_value(cb_data, x)` is integer (to some tolerance).
        #return
    else
        @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        # `callback_value(cb_data, x)` might be fractional or integer.
    end

        x_val = zeros((2))
        for i in 1:2
            x_val[i] = callback_value(cb_data, x[i])
        end

        y_val = callback_value(cb_data, y)
        if sum(x_val)+y_val > 13 + 1e-6
            con = @build_constraint(sum(x) +y <= 13)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    
end
MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)


JuMP.optimize!(model)
