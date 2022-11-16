using Catlab, Catlab.Graphs, Catlab.Graphics, Catlab.CategoricalAlgebra
using Catlab.Theories, Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = AlgebraicRewriting.homomorphism

function get_valid_op1s(deca_source, varID)
    # skip_ops = Set([:∂ₜ])
    indices = incident(deca_source, varID, :tgt)
    return filter!(x -> deca_source[x, :op1] != :∂ₜ, indices)
end
  
function get_target_indices(deca_source)
    targetVars = []
    for var in parts(deca_source, :Var)
    op1Count = length(get_valid_op1s(deca_source, var))
    op2Count = length(incident(deca_source, var, :res))
    sumCount = length(incident(deca_source, var, :sum))

    tot = op1Count + op2Count + sumCount
    if(tot >= 2)
        append!(targetVars, var)
    end
    end

    return targetVars
end
  
function get_preprocess_indices(deca_source)
    targetOp2 = []
    targetSum = []

    targetVars = get_target_indices(deca_source)

    for var in targetVars
    append!(targetOp2, incident(deca_source, var, :res))
    append!(targetSum, incident(deca_source, var, :sum))
    end

    return targetOp2, targetSum
end
  
function preprocess_rewrite(deca_source)
    targetOp2, targetSum = get_preprocess_indices(deca_source)

    # If we don't need to preprocess then don't
    if(length(targetOp2) == 0 && length(targetSum) == 0)
    return deca_source
    end

    LHS = []
    RHS = []

    SuperMatch = []
    SuperVarMap = Vector{Int64}()
    SuperOp2Map = Vector{Int64}()
    SuperSigmaMap = Vector{Int64}()
    SuperSummandMap = Vector{Int64}()

    serial = 0
    # Process all of the target rewrites for op2
    for opID in targetOp2

    vars = [deca_source[opID, :proj1], deca_source[opID, :proj2], deca_source[opID, :res]]
    types = deca_source[vars, :type]
    names = deca_source[vars, :name]
    op2name = deca_source[opID, :op2]

    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = 3
        type = types
        name = names

        Op2 = 1
        proj1 = [1]
        proj2 = [2]
        res = [3]
        op2 = op2name
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = 3
        type = types
        name = names
    end

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = 4
        type = vcat(types, types[end])
        name = vcat(names, Symbol("Temp_", serial))

        Op1 = 1
        src = [4]
        tgt = [3]
        op1 = [:temp]

        Op2 = 1
        proj1 = [1]
        proj2 = [2]
        res = [4]
        op2 = op2name
    end

    serial += 1

    L = ACSetTransformation(I, Match, Var = 1:3)
    R = ACSetTransformation(I, Sub, Var = 1:3)

    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    push!(SuperOp2Map, opID)
    end
    
    # Process all of the target rewrites for sums
    for sumID in targetSum
    summandIDs = incident(deca_source, sumID, :summation)
    vars = vcat(deca_source[summandIDs, :summand], deca_source[sumID, :sum])
    types = deca_source[vars, :type]
    names = deca_source[vars, :name]
    
    rewrite_size = length(vars)
    nary = rewrite_size - 1

    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = rewrite_size
        type = types
        name = names

        Σ = 1
        sum = [rewrite_size]

        Summand = nary
        summand = 1:nary
        summation = fill(1, nary)
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = rewrite_size
        type = types
        name = names
    end

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = rewrite_size + 1
        type = vcat(types, types[end])
        name = vcat(names, Symbol("Temp_", serial))

        Op1 = 1
        src = [rewrite_size + 1]
        tgt = [rewrite_size]
        op1 = [:temp]

        Σ = 1
        sum = [rewrite_size + 1]

        Summand = nary
        summand = 1:nary
        summation = fill(1, nary)
    end

    serial += 1

    L = ACSetTransformation(I, Match, Var = 1:rewrite_size)
    R = ACSetTransformation(I, Sub, Var = 1:rewrite_size)

    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    push!(SuperSigmaMap, sumID)
    append!(SuperSummandMap, summandIDs)
    end

    # Combine all rules in parallel and apply
    rule = Rule(oplus(LHS), oplus(RHS))
    m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op2 = SuperOp2Map, Σ = SuperSigmaMap, Summand = SuperSummandMap)

    rewrite_match(rule, m)
end
  
function rewrite_decapode(deca_source)
    # Just for now, I'll get it working with op1 only
    # Considering op2 and summations will make this significantly more difficult

    targetVars = get_target_indices(deca_source)

    if(length(targetVars) == 0)
    return deca_source
    end

    LHS = []
    RHS = []

    SuperMatch = []
    SuperVarMap = Vector{Int64}()
    SuperOp1Map = Vector{Int64}()

    varSerial = 0
    sumSerial = 0
    for varID in targetVars
    targetOp1 = get_valid_op1s(deca_source, varID)
    vars = vcat(deca_source[targetOp1, :src], varID)

    num_nodes_match = length(vars)
    nary_of_rewrite = num_nodes_match - 1

    result_index = num_nodes_match
    sum_index = 2 * result_index

    variable_types = deca_source[vars, :type]
    variable_var =  deca_source[vars, :name]
    variable_op1 = deca_source[targetOp1, :op1]

    # Cannot work with undefined types
    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = num_nodes_match
        type = variable_types
        name = variable_var

        # This will probably break for rewrites including 
        # Non-Op1 rewrites
        Op1 = nary_of_rewrite
        src = 1:nary_of_rewrite
        tgt = fill(result_index, nary_of_rewrite)
        op1 = variable_op1
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = num_nodes_match
        type = variable_types
        name = variable_var
    end 

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
        Var = 2 * num_nodes_match
        type = vcat(variable_types, variable_types)
        name = vcat(variable_var, map(x -> Symbol("••", varSerial + x), 1:nary_of_rewrite), [Symbol("••sum", sumSerial)])
        Op1 = nary_of_rewrite + 1
        src = vcat(1:nary_of_rewrite, sum_index)
        tgt = vcat(num_nodes_match+1:sum_index-1, [result_index])
        op1 = vcat(variable_op1, Symbol(:avg, nary_of_rewrite))
        Σ = 1
        sum = [sum_index]
        Summand = nary_of_rewrite
        summand = num_nodes_match+1:sum_index-1
        summation = fill(1, nary_of_rewrite)
    end

    varSerial += nary_of_rewrite
    sumSerial += 1

    L = ACSetTransformation(I, Match, Var = 1:num_nodes_match);
    R = ACSetTransformation(I, Sub, Var = 1:num_nodes_match);
    
    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    append!(SuperOp1Map, targetOp1)
    end

    rule = Rule(oplus(LHS), oplus(RHS))

    m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op1 = SuperOp1Map)
    rewrite_match(rule, m)
end
  
function average_rewrite(deca_source::SummationDecapode)
    return rewrite_decapode(preprocess_rewrite(deca_source))
end