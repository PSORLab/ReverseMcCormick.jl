using PkgBenchmark
results = benchmarkpkg("ReverseMcCormick")
show(results)

#=
# specify tag and uncommit to benchmark versus prior tagged version
tag =
results = judge("ReverseMcCormick", tag)
show(results)
=#

export_markdown("results.md", results)
