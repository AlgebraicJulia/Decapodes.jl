using Aqua, Decapodes
Aqua.test_all(
  Decapodes, ambiguities=false,
  deps_compat=(ignore=[:Markdown, :Random, :Test],),
)
