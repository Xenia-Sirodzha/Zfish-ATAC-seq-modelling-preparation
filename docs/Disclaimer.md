The output of the Rmd is not sufficiently rendered with the following settings in the R notebook:

```
---
title: "Title Name"
date: "2025-02-06"
author: Xenia Sirodzha
output:
  github_document: 
    toc: yes
    toc_depth: 4
---
```
 While an .md file is rendered, graphs are not sufficiently rendered. Therefore, the code is only available as a .Rmd file, and the rendered output is only available as a downloadable .html file. These files are in the folders `Rmds` and `htmls` here in `docs`.
There have been no successful attempts to fix this issue yet. 
