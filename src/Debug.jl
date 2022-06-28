module Debug
  using CombinatorialSpaces
  using CombinatorialSpaces: volume
  using Catlab.Graphics
  using Catlab.WiringDiagrams
  using Catlab.CategoricalAlgebra
  using Catlab.Theories
  using GeometryBasics
  using ...CairoMakie
  export sim_key, get_wire, draw_wire

  function sim_key(dwd; kw...)
    n_dwd = deepcopy(dwd)
    outer_ports = nparts(n_dwd.diagram, :OuterInPort)
    n_dwd.diagram[:outer_in_port_type] .= "inp"
    n_dwd.diagram[:out_port_type] .= (1:nparts(n_dwd.diagram, :OutPort))
    to_graphviz(n_dwd; labels=true, kw...)
  end

  function get_wire(dwd, vf, u0, wire; p=[], t=0)
    vf(copy(u0), u0, p, t)
    return vf.mem.contents[wire]
  end

  function draw_wire(s, sd, dwd, vf, u0, wire; p=[], t=0, colorrange=nothing, xlim=nothing, ylim=nothing, axisaspect=nothing, kw...)
    w_data = get_wire(dwd, vf, u0, wire; p=p, t=t)
    w_type = supertype(typeof(dwd.diagram[wire, :out_port_type][:type]))
    plot(w_type, s, sd, w_data, colorrange, xlim, ylim, axisaspect; kw...)
  end

  function plot(::Type{ObExpr{:Form0}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
    fig = Figure()
    ax, ob = mesh(fig[1,1], s, color=vals)
    if !isnothing(colorrange)
      ob.colorrange = colorrange
    end
    if !isnothing(xlim)
      xlims!(ax, xlim)
    end
    if !isnothing(ylim)
      ylims!(ax, ylim)
    end
    if !isnothing(axisaspect)
      ax.aspect = AxisAspect(axisaspect)
    end
    Colorbar(fig[1, 2], ob)
    fig, ax, ob
  end

  function plot(::Type{ObExpr{:DualForm2}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
    vals ./= [sum(dual_volume(Val{2}, sd, elementary_duals(Val{0},sd,v))) for v in 1:nv(s)]
    plot(ObExpr{:Form0}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
  end

  function plot(::Type{ObExpr{:DualForm1}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
    vals .*= [volume(Val{1},sd,e) / sum(dual_volume(Val{1}, sd, elementary_duals(Val{1},sd,e)))
              for e in 1:ne(s)]
    plot(ObExpr{:Form1}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
  end

  function plot(::Type{ObExpr{:DualForm0}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
    vals .*= [volume(Val{2},sd,t) for t in 1:ntriangles(s)]
    plot(ObExpr{:Form2}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
  end

  function plot(::Type{ObExpr{:Form1}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; use_arrows=false, n_arrows = 100, kw...)
    vals ./= [volume(Val{1},sd,e) for e in 1:ne(s)]
    if isnothing(colorrange)
      colorrange = (0.0, maximum(abs.(vals)))
    end
    fig = Figure()
    ax, ob = if use_arrows
      orient_vals = vals .* [eo ? 1 : -1 for eo in s[:edge_orientation]]
      signs = sign.(orient_vals)
      locs = s[s[:∂v1], :point] .* (0.5 .- (signs * 0.35)) .+ s[s[:∂v0], :point] .* (0.5 .+ (signs * 0.35))
      mag = ((s[s[:∂v1], :point] .- s[s[:∂v0], :point]) * 0.5) .* signs
      inds = collect(1:ne(s))
      if !isnothing(xlim)
        filter!(i -> xlim[1] <= s[s[i, :∂v0], :point][1] <= xlim[2], inds)
      end
      if !isnothing(ylim)
        filter!(i -> ylim[1] <= s[s[i, :∂v0], :point][2] <= ylim[2], inds)
      end
      inds = inds[1:min(n_arrows, length(inds))]

      arrows(fig[1,1], locs[inds], mag[inds], color=abs.(vals[inds]), colorrange=colorrange)
    else
      linesegments(fig[1,1], s, color=abs.(vals), colorrange=colorrange)
    end
    if !isnothing(xlim)
      xlims!(ax, xlim)
    end
    if !isnothing(ylim)
      ylims!(ax, ylim)
    end
    if !isnothing(axisaspect)
      ax.aspect = AxisAspect(axisaspect)
    end
    Colorbar(fig[1, 2], limits=colorrange)
    fig, ax, ob
  end

  function plot(::Type{ObExpr{:Form2}}, s, sd, vals, colorrange, xlim, ylim, axisaspect; kw...)
    # Split mesh into component triangles
    m = GeometryBasics.Mesh(s)
    x = faces(m)
    m_points = m.position[vcat([[t[1],t[2],t[3]] for t in x]...)]
    m_faces = TriangleFace{Int}[[((t-1) * 3) .+ (1,2,3) for t in  1:length(x)]...]
    new_m = GeometryBasics.Mesh(Point{3, Float64}[m_points...], m_faces)
    fig = Figure()
    ax, ob = mesh(fig[1,1], new_m, color=vcat([[v,v,v] for v in vals]...))
    if !isnothing(colorrange)
      ob.colorrange = colorrange
    end
    if !isnothing(xlim)
      xlims!(ax, xlim)
    end
    if !isnothing(ylim)
      ylims!(ax, ylim)
    end
    if !isnothing(axisaspect)
      ax.aspect = AxisAspect(axisaspect)
    end
    Colorbar(fig[1, 2], ob)
    fig, ax, ob
  end
end
