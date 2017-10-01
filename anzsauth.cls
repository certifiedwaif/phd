%
%---------------------------------------------------------------------------
% This file is a slight modification of an earlier (2013) version.
% These modifications were made by Rolf Turner (Technical Editor
% of ANZJS) as of 19/05/2016.  The modifications are designed 
% to accommodate a few changes that have been made to the Journal's
% policies and protocols in respect of presentation style.
%
% This class file is made available for use by authors who wish to
% prepare an article for publication in the Australian & New Zealand
% Journal of Statistics published by Wiley Publishing Asia Pty Ltd.
%
%This class file is provided on an `as is'  basis, without warranties of any
%kind, either expressed or implied, including but not limited to warranties of
%title, or  implied warranties of merchantablility or fitness for  a
%particular purpose. There will be no duty on the author[s] of the software
%or Wiley Publishing Asia Pty Ltd. to correct any errors or defects in the software.
%Any statutory rights you may have remain unaffected by your acceptance of
%these rules of use.
%---------------------------------------------------------------------------
%
% Origonal version created by Alistair Smith, Sunrise Setting Ltd, 19 April 2013
%
% anzsauth.cls --- For Aust. N. Z. J. Stat.
% 6/1/14 change to documentation
%
%\def\update{2014/01/06 v1.01}
% Changes made by Rolf Turner 19/05/2016.
\def\update{2016/19/05 v1.02}

\newcommand{\journalname}{AUSTRALIAN \& NEW ZEALAND JOURNAL OF STATISTICS}
\newcommand{\journalnamelc}{Australian \& New Zealand Journal of Statistics}
\newcommand{\journalabb}{Aust. N.~Z. J. Stat.}
\newcommand{\journalclass}{anzsauth.cls}
\newcommand{\journalclassshort}{anzsauth}
\newcommand{\DOI}{XXX}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{anzsauth}[\update\ \journalclass]

%\newcommand\hmmax{0}

\newif\if@timesfont
\DeclareOption{times}{%
\@timesfonttrue}

\newif\if@doublespace
\DeclareOption{doublespace}{%
\@doublespacetrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
\if@timesfont
\RequirePackage{times}
\fi
\if@doublespace
\RequirePackage[onehalfspacing]{setspace}
\fi

\RequirePackage{graphicx}
\RequirePackage{pifont,latexsym,ifthen,rotating,calc,textcase,booktabs,color}
\RequirePackage{amsfonts,amssymb,amsbsy,amsmath,amsthm}
%\RequirePackage{bm}
\RequirePackage[errorshow]{tracefnt}

\@twosidetrue
\flushbottom
\frenchspacing

\textwidth 130mm
\textheight 573pt%
%\setlength\columnsep{24pt}

%Trim sizes
\setlength\voffset{-1in}
\setlength\hoffset{-1in}
\topmargin -1mm
\setlength\oddsidemargin{20.7mm}%back margin on odd pages
\setlength\evensidemargin{20.7mm}%fore margin on even pages
\setlength{\pdfpagewidth}{171.4mm}%
\setlength{\pdfpageheight}{245.9mm}%
\setlength\paperwidth{171.4mm}
\setlength\paperheight{245.9mm}
%Needed to set PDF page size
\special{papersize=171.4mm,245.9mm}

\parskip \z@
\parindent 2em
\headheight 50pt
\headsep    14pt 
\footskip 24pt 

\hyphenpenalty=1000
\pretolerance=8000
\tolerance=9500
\hbadness=8000
\vbadness=9000
\displaywidowpenalty=0
\clubpenalty=10000
\widowpenalty=10000
\lefthyphenmin=3%
\righthyphenmin=3%
\brokenpenalty=10000%

\thinmuskip = 3mu
\medmuskip = 4mu
\thickmuskip = 5mu

\setcounter{topnumber}{10}
\def\topfraction{1}
\setcounter{bottomnumber}{10}
\def\bottomfraction{0.8}
\setcounter{totalnumber}{10}
\def\textfraction{0}
\renewcommand{\floatpagefraction}{0.95}
\setcounter{dbltopnumber}{10}
\renewcommand{\dblfloatpagefraction}{0.95}
\renewcommand{\dbltopfraction}{1}

\renewcommand{\normalsize}{\fontsize{10.0}{12.5pt}\selectfont}
\renewcommand{\small}{\fontsize{9}{10pt}\selectfont}
\renewcommand{\footnotesize}{\fontsize{8.5}{9pt}\selectfont}
\renewcommand{\scriptsize}{\fontsize{7.5}{8.5pt}\selectfont}
\renewcommand{\tiny}{\fontsize{6.5}{7pt}\selectfont}
\renewcommand{\large}{\fontsize{11}{12pt}\selectfont}
\renewcommand{\Large}{\fontsize{14}{18pt}\selectfont}
\renewcommand{\LARGE}{\fontsize{17}{22pt}\selectfont}
\renewcommand{\huge}{\fontsize{20}{25pt}\selectfont}
\renewcommand{\Huge}{\fontsize{25}{30pt}\selectfont}

\newcommand{\titlesize}{\fontsize{15.0}{16pt}\selectfont}
\newcommand{\tabsize}{\fontsize{9}{9.5pt}\selectfont}
\newcommand{\headersize}{\fontsize{8.0}{9pt}\selectfont}

\newbox\absbox
\def\abstract{\lrbox\absbox\minipage{116mm}%
  \normalfont\small%
  \section*{Summary}%
  }
\def\endabstract{\endminipage\endlrbox}

\def\keywords#1{%
  \gdef\@keywords{\begin{minipage}{116mm}{\footnotesize\sl Key words:}\ \parbox[t]{98mm}{#1}\end{minipage}}}
\let\@keywords\@empty

\skip\footins 20pt plus  8pt
%\gdef\footnoterule{}
\def\footnoterule{\kern-3\p@
  \hrule \@width 144pt \kern 2.6\p@}

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\long\def\@makefntext#1{\parindent 1em%
\noindent{$\m@th^{\@thefnmark}$}#1}

\def\corrauth{\footnotemark[1]\gdef\@corraddr{%
 \footnotetext{\footnotemark[1]Author to whom correspondence should be addressed.\stepcounter{footnote}}}}
\let\@corraddr\@empty

\def\affiliation#1{%
  \gdef\@affiliation{#1}}
\let\@affiliation\@empty

\def\address#1{%
  \gdef\@address{%
    \footnotetext[0]{#1}}}
\let\@address\@empty

\def\acks#1{%
  \gdef\@acks{%
   {\footnotetext[0]{{\sl Acknowledgments.} #1}}}}
\let\@acks\@empty

\def\ack#1{%
  \gdef\@ack{%
   {\footnotetext[0]{{\sl Acknowledgment.} #1}}}}
\let\@ack\@empty

\def\addressnum#1{${}^{#1}$}
%\def\affil#1{${}^{#1}$}
\def\comma{${}^{\text{,}}$}

\renewcommand\maketitle{\par
  \begingroup
     \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
       \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{title}\label{FirstPage}\@corraddr\@address\@acks\@ack%\@cgs\@cgsn
  \endgroup
  %\setcounter{footnote}{0}%
  \global\let\affiliation\relax
  \global\let\address\relax
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\@affiliation\@empty
  \global\let\@address\@empty
  %\global\let\corraddr\relax
  \global\let\acks\relax
  \global\let\ack\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{\vspace*{-18pt}%
\null%
\begin{center}
    {\normalsize{\bf\@title} \par}%
    \vskip 1.5em %
    {\large
      \lineskip .5em%
      \textsc{\@author}
      \par}%
      \vskip 7pt
{\large
      \lineskip .5em%
      %
        {\sl\@affiliation}
      \par}%
    \vskip 19pt%
    {\noindent\usebox\absbox\par}
    {\vspace{5pt}%
      %
      {\noindent\footnotesize\@keywords}\par}
      \end{center}
      \vspace{-10pt}
        \par%
  }
 
\def\startpage{\pageref{FirstPage}}
\def\endpage{\pageref{LastPage}}
\def\volumeyear{0000}
\def\volumenumber{00}

\gdef\copyrightline{\copyright\ \volumeyear\ Australian Statistical Publishing Association Inc.}
\def\runningheads#1#2{\markboth{\uppercase{#1}}{\uppercase{#2}}}

\def\ps@title{%
\def\@oddhead{%
\parbox[t]{\textwidth}{%
{\footnotesize{\sl\journalabb}\ \volumeyear\hfill doi: 10.1111/j.1467-842X.\DOI}}}
\let\@evenhead\@oddhead
\def\@oddfoot{\parbox[t]{\textwidth}{%
{\scriptsize\copyrightline\ Published by Wiley Publishing Asia Pty Ltd.}\\
{\scriptsize\textsl{Prepared using \textsf{\journalclass} [Version: \update]}}}}
\let\@evenfoot\@oddfoot}

\def\ps@wpage{
\let\@mkboth\@gobbletwo
\def\@evenhead{\footnotesize \thepage\hfill\footnotesize{\leftmark}\hfill\footnotesize\phantom{\thepage}}
\def\@oddhead{\footnotesize \phantom{\thepage}\hfill\footnotesize{\rightmark}\hfill\footnotesize\thepage}
\def\@evenfoot{\parbox[t]{\textwidth}{{\scriptsize \copyrightline}\\ %
\scriptsize\textsl{Prepared using \textsf{\journalclass}}}}
\def\@oddfoot{\@evenfoot}
}

\renewcommand{\@seccntformat}[1]{{\csname the#1\endcsname.}\hspace{0.5em}}

\newdimen\@bls                              
\@bls=\baselineskip

\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {1.3\@bls  plus .3\@bls minus .1\@bls}%
                                   {8pt}%\@afterindentfalse}%
                                   {\centering\normalfont\normalsize\bf}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {\@bls plus .3\@bls minus .1\@bls}%
                                     {6pt}%\@afterindentfalse}%
                                     {\normalfont\normalsize\bf\raggedright}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {\@bls plus .3\@bls minus .1\@bls}%
                                     {6pt}%\@afterindentfalse}%
                                     {\normalfont\normalsize\bfseries\slshape}}

\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\itemsep 0pt\parsep 0pt
         \def\makelabel##1{\hss\llap{##1}}}\fi}

\let\endenumerate =\endlist

\def\itemize{\ifnum \@itemdepth >3 \@toodeep\else \advance\@itemdepth \@ne
\edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
\list{\csname\@itemitem\endcsname}{\itemsep 0pt\parsep 0pt
\def\makelabel##1{\hss\llap{##1}}}\fi}

\let\enditemize =\endlist

%\renewcommand{\thetable}{\Roman{table}}

\usepackage{caption}
\DeclareCaptionLabelSeparator{jwperiod}{.\hspace*{1.0ex}}
%\DeclareCaptionFont{capsize}{\tabsize}
\DeclareCaptionFormat{anzsformat}{\centerline{#1#2}\\[3pt]#3}
\captionsetup[figure]{font=small,labelfont=rm,labelsep=jwperiod,
   justification=justified,singlelinecheck=true}
\captionsetup[table]{format=anzsformat,position=top,labelsep=none,
   font=small,labelfont={sc,large},textfont=sl,justification=justified}
\DeclareCaptionStyle{oneline}{format=anzsformat,position=top,labelsep=none,
   font=small,labelfont={sc,large},textfont=sl,justification=centering}

\def\@begintheorem#1#2[#3]{%
  \deferred@thm@head{\the\thm@headfont \thm@indent
    \@ifempty{#1}{\let\thmname\@gobble}{\let\thmname\@iden}%
    \@ifempty{#2}{\let\thmnumber\@gobble}{\let\thmnumber\@iden}%
    \@ifempty{#3}{\let\thmnote\@gobble}{\let\thmnote\@iden}%
    \thm@swap\swappedhead\thmhead{#1}{#2\the\thm@headpunct}{#3}%
    %\the\thm@headpunct
    \thmheadnl % possibly a newline.
    \hskip\thm@headsep
  }%
  \ignorespaces}

\def\thmhead@plain#1#2#3{%
  \thmname{#1}\thmnumber{\@ifnotempty{#1}{ }{#2}}%
  \thmnote{ {\the\thm@notefont #3. }}}

\newtheoremstyle{wiley}
{6pt plus 2pt minus 2pt}% space above
{6pt plus 2pt minus 2pt}% space below
{\sl}% Body font
{}% Indent amount
{\bfseries}% Theorem head font
{.}% Punctuation after theorem head
{0.75em}% Space after theorem head
{}% Theorem head spec

\theoremstyle{wiley}

\renewenvironment{proof}[1][\proofname]{\par
  %AS\pushQED{\qed}%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        \bfseries
    #1\@addpunct{.}]\ignorespaces
}{%
  %AS\popQED
  \endtrivlist\@endpefalse
}

\renewenvironment{thebibliography}{%
\section*{\refname}
\list{}{\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=1\p@\labelsep=0\p@\itemindent=-20\p@
\labelwidth=0\p@\leftmargin=20\p@
}\refsize\rm
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}{\endlist}

\usepackage{natbib}
\bibpunct{(}{)}{;}{}{}{,}
\setlength{\bibsep}{1pt}
\setlength{\bibhang}{20pt}
\renewcommand{\bibfont}{\footnotesize}
\renewcommand\refname{{\normalfont\sl References}}

\AtEndDocument{%
\label{LastPage}}

\pagestyle{wpage}
\normalsize
\sloppy
