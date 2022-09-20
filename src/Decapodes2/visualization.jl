### Drawing of Decapodes
import Catlab.Graphics.Graphviz

#This could probably be made neater
function Catlab.Graphics.to_graphviz(d::NamedDecapode)::Graphviz.Graph
    #Similar to the to_graphviz in other implementations
    gv_name(v::Int) = "n$v"
    
    gv_path(e::Int) = [gv_name(d[:src][e]), gv_name(d[:tgt][e])]

    gp_name(p::Int) = "p$p"
    gp_proj1(p::Int) = [gp_name(p), gv_name(d[:proj1][p])]
    gp_proj2(p::Int) = [gp_name(p), gv_name(d[:proj2][p])]
    gp_projRes(p::Int) = [gp_name(p), gv_name(d[:res][p])]

    stmts = Graphviz.Statement[]

    reg_to_sub = Dict('0'=>'₀', '1'=>"₁", '2'=>'₂', '3'=>'₃', '4'=>'₄',
    '5'=>'₅', '6'=>'₆','7'=>'₇', '8'=>'₈', '9'=>'₉')

    toSub(digit::Char) = haskey(reg_to_sub, digit) ? reg_to_sub[digit] : digit

    #For variables, label grabs the stored variable name and its type and concatenate
    #label assumes dimension is single digit
    for v in parts(d, :Var)
        vertex_name = String(d[:name][v]) * ":Ω" * toSub(last(String(d[:type][v])))
        push!(stmts, Graphviz.Node(gv_name(v), Dict(:label=>vertex_name)))
    end

    #For unary ops, label mashes together all func symbol names into one string
    for e in parts(d, :Op1)
        #add composition symbols?
        edge_name = join(String.(d[:op1][e]))
        push!(stmts, Graphviz.Edge(gv_path(e), Dict(:label=>edge_name)))
    end

    #For binary ops, make temp product object, drop projections and drop result with op name
    for p in parts(d, :Op2)
        proj_space_name = "Ω" * toSub(last(String(d[:type][d[:proj1][p]]))) * "×" * "Ω" * toSub(last(String(d[:type][d[:proj2][p]])))
        push!(stmts, Graphviz.Node(gp_name(p), Dict(:label=>proj_space_name)))

        push!(stmts, Graphviz.Edge(gp_proj1(p), Dict(:label=>"proj₁", :style=>"dashed")))
        push!(stmts, Graphviz.Edge(gp_proj2(p), Dict(:label=>"proj₂", :style=>"dashed")))

        res_name = String(d[:op2][p])
        push!(stmts, Graphviz.Edge(gp_projRes(p), Dict(:label=>res_name)))
    end

    #Need to add user access for more customizability later
    Graphviz.Graph("G", true, "neato", stmts, Dict(), Dict(:shape=>"oval"), Dict())
end
