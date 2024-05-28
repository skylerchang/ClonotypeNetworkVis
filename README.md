# Immune Repertoire Clonotype Network Dashboard

## Introduction
This dashboard is designed to provide a clean, simple and accessible interface to explore and virsualize highly similar clonotypes. 

Currently, the dashboard includes the following tabs and features:





<video width="320" height="240" controls>
<source src="https://github.com/skylerchang/ClonotypeNetworkVis/main/video_tutorials/Part1a.MP4" type="video/MP4">
</video>


This dashboard was built using many great tools in the R ecosystem. Thanks to all of the developers of these open source packages:

- [shiny]
- [rtweet]
- [shinydashboard]
- [plotly]
- [tidyverse]
- [shinycssloaders]
- [DT]

...and many more. For a full list of project dependencies, see [deps.yaml](deps.yaml).

I also built a few things to make this work, including:

- [gathertweet] - A command line tool for gathering tweets from Twitter search streams. Removes the boilerplate of collecting Twitter data and plays nicely with `cron`.

- [shinyThings] - [shiny] modules for pagination and dropdown buttons.

- [shinydashboard][shinydashboard-fork] (fork) - I forked [shinydashboard] to add a few features I needed to make this work the way I wanted.

- [adminlte-ocean-next] - An [AdminLTE] dashboard color theme.
    
---

This dashboard was built by [Garrick Aden-Buie][garrick-home] and is released under an [MIT license][mit-license].

You are welcome to re-use and/or customize this dashboard! If you do, I kindly request that you provide credit and link back to the [source repo][repo] or my [personal webpage][garrick-home]. Thank you!
