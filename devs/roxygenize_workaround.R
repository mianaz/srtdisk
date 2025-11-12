#!/usr/bin/env Rscript
# Workaround for roxygen2 parser_setGeneric bug
# This patches the parser to skip setGeneric processing but preserve metadata

# Monkey-patch roxygen2's parser_setGeneric function
assignInNamespace(
  "parser_setGeneric",
  function(call, env, block) {
    # Extract the generic name from the call
    name <- as.character(call$name)

    # Create a minimal object that roxygen2 can process
    # We create a simple function object instead of trying to get the S4 generic
    value <- function() {}

    # Return an object with the name
    roxygen2:::object(value, name, "function")
  },
  ns = "roxygen2"
)

# Now run roxygenize
message("Running roxygenize with patched parser...")
roxygen2::roxygenize()
message("Documentation generation complete!")
