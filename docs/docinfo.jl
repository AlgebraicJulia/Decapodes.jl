module DocInfo

using Dates

mutable struct Info
  start_time::DateTime
  finish_time::DateTime
end

function Info()
  Info(now(), now())
end

function get_elapsed(info::Info)
  info.finish_time - info.start_time
end

function get_report(info::Info)
  info.finish_time = now()
  elapsed = get_elapsed(info)
  elapsed_sec = round(elapsed, Dates.Second(1))
  @info "Page built in $(elapsed_sec)."
  @info "This page was last built at $(info.finish_time)."
end

draw(deca) = to_graphviz(deca, box_labels=:name, junction_labels=:variable, prog="circo")

end
