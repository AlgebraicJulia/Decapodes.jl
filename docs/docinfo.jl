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
  @info "Documentation built in $(elapsed_sec)."
  @info "Documentation built at $(info.finish_time)."
end

end