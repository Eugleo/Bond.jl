using Documenter
using Bond

makedocs(
    sitename = "Bond",
    format = Documenter.HTML(),
    modules = [Bond]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
