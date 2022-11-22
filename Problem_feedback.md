## Questions about RankCompV3 feedback 

**The following is a summary of the current user feedback problems and solutions, I hope to help you.**

##### 1. ERROR: UndefVarError: writeshortest is not defined.

This is due to Parsers 2.5.0's writeshortest function being unavailable in early November. You can go to add Parsers@2.4.2 for version rollback compatibility.

[Playing with ideas on how to better componentize internal parsing by quinnj · Pull Request #127 · JuliaData/Parsers.jl (github.com)](https://github.com/JuliaData/Parsers.jl/pull/127)
