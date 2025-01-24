using PicardLefschetz
using Documenter

DocMeta.setdocmeta!(PicardLefschetz, :DocTestSetup, :(using PicardLefschetz); recursive=true)

makedocs(;
    modules=[PicardLefschetz],
    authors="Anne Weber <anne.weber@mailbox.org> and contributors",
    sitename="PicardLefschetz.jl",
    format=Documenter.HTML(;
        canonical="https://anneaux.github.io/PicardLefschetz.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" =>"examples.md"
    ],
)

deploydocs(;
    repo="github.com/anneaux/PicardLefschetz.jl",
    devbranch="main",
)
