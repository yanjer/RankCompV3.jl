## Questions about RankCompV3 feedback 

**The following is a summary of the current user feedback problems and solutions, I hope to help you.**

##### 1. ERROR: UndefVarError: writeshortest is not defined.

This is because the writeshortest function was disabled in Parsers 2.5.0 in early November. You can use  Parsers2.4.2 for version back compatibility.

```julia
add Parsers@2.4.2
```

[Playing with ideas on how to better componentize internal parsing by quinnj · Pull Request #127 · JuliaData/Parsers.jl (github.com)](https://github.com/JuliaData/Parsers.jl/pull/127)
