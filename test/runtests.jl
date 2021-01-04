using Test
using DiffinDiffsBase

include("testutils.jl")

const tests = [
    "utils",
    "terms",
    "treatments",
    "parallels"
]

printstyled("Running tests:\n", color=:blue)

for test in tests
    include("$test.jl")
    println("\033[1m\033[32mPASSED\033[0m: $(test)")
end
