\documentclass[10pt, usenames,dvipsnames, xcolor=table]{beamer}
\usepackage[utf8]{inputenc}
\usepackage{utopia}
\usepackage{lmodern}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{nicefrac}
\usepackage{parskip}
\usepackage{media9}
\definecolor{UniBlue}{RGB}{0, 81, 149}
\setbeamercolor{structure}{fg=UniBlue}
\usetheme{Darmstadt}
\setbeamersize{text margin left=3mm,text margin right=3mm} 
\usepackage{gensymb}
\usepackage{mathtools}
\usepackage{hyperref}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{cancel}
\usepackage{multimedia}
\usepackage{resizegather}
\setcounter{MaxMatrixCols}{20}
\graphicspath{{Figures/}}

%------------------------------------------------------------
%This block of code defines the information to appear in the
%Title page

\makeatother
\setbeamertemplate{footline}
{
  \leavevmode%
  \hbox{%
  \begin{beamercolorbox}[wd=.2\paperwidth,ht=2.25ex,dp=1ex,center]{author in head/foot}%
    \usebeamerfont{author in head/foot}\insertshortauthor
  \end{beamercolorbox}%
  \begin{beamercolorbox}[wd=.8\paperwidth,ht=2.25ex,dp=1ex,center]{title in head/foot}%
    \usebeamerfont{title in head/foot}\insertshorttitle\hspace*{3em}
    \insertframenumber{} / \inserttotalframenumber\hspace*{1ex}
  \end{beamercolorbox}}%
  \vskip0pt%
}
\makeatletter
\setbeamertemplate{navigation symbols}{}
