% LaTeX2e thesis style for Maths at the University of Sydney.
% It sets up a title page, loads up the amsbook class file and
% then redefines a number of things which I don't like in amsbook.

% Andrew Mathas May, 2002

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{usydthesis}
  [2002/05/10 University of Sydney thesis class file - Andrew Mathas]

% options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\SMSsection{Pure Mathematics}
\DeclareOption{pure}{\def\SMSsection{Pure Mathematics}}
\DeclareOption{applied}{\def\SMSsection{Applied Mathematics}}
\DeclareOption{stats}{\def\SMSsection{Statistics}}

\newif\ifHandBook\HandBookfalse
\DeclareOption{handbook}{\HandBooktrue}

\def\MasterClass{amsbook}

% degree types
\def\@thesistype{An essay}
\def\@degree{something cool}
\DeclareOption{BA}{\def\@thesistype{An essay}
                   \def\@degree{B.A.~(Honours)}
}
\DeclareOption{BSc}{\def\@thesistype{An essay}
                    \def\@degree{B.Sc.~(Honours)}
}
\DeclareOption{MA}{\def\@thesistype{A thesis}
                    \def\@degree{Master of Arts}
}
\DeclareOption{MSc}{\def\@thesistype{A thesis}
                    \def\@degree{Master of Science}
}
\DeclareOption{PhD}{\def\@thesistype{A thesis}
                    \def\@degree{Doctor of Philosophy}
}
\newif\ifBook\Bookfalse
\DeclareOption{book}{\Booktrue\def\MasterClass{book}}

\newif\ifSMStimes\SMStimestrue
\DeclareOption{notimes}{\SMStimesfalse}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}\MasterClass}
\ProcessOptions\relax

\LoadClass[12pt]\MasterClass
\usepackage{amsmath,amsthm}

\hfuzz=4pt
\vfuzz=4pt

% packages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the University crest
\usepackage{graphicx}

% amsmath
\usepackage{amssymb,amsmath}

% times fonts look better
\ifSMStimes
  \usepackage{times}
  \usepackage[T1]{fontenc}
  \usepackage{textcomp}
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we now change a few things in the standard amsbook setup

% sections are numbered within the chapters
\numberwithin{section}{chapter}

% chapters don't need to start on the left hand page
\def\chapter{%
  \clearpage
  \thispagestyle{plain}\global\@topnum\z@
  \@afterindenttrue \secdef\@chapter\@schapter}

% section numbers in appendices shouldn't have a dot in them
\def\appendix{\par
  \c@chapter\z@ \c@section\z@
  \let\chaptername\appendixname
  \def\thechapter{\@Alph\c@chapter}
  \renewcommand\thesection       {\thechapter\arabic{section}}
}

\ifBook\relax
\else
 %%put dotfill in toc entries and make chapter entries bold 
  \def\contentsname{\normalfont\scshape{Contents}}
  \def\@tocline#1#2#3#4#5#6#7{\relax
  \ifnum #1>\c@tocdepth % then omit
  \else
    \par \addpenalty\@secpenalty\addvspace{#2}%
    \begingroup \hyphenpenalty\@M
    \@ifempty{#4}{%
      \@tempdima\csname r@tocindent\number#1\endcsname\relax
    }{%
      \@tempdima#4\relax
    }%
    \parindent\z@ \leftskip#3\relax \advance\leftskip\@tempdima\relax
    \rightskip\@pnumwidth plus4em \parfillskip-\@pnumwidth
    {#5\leavevmode\hskip-\@tempdima #6}\nobreak\relax
    \dotfill\hbox to\@pnumwidth{\@tocpagenum{#5{#7}}}\par
    \nobreak
    \endgroup
  \fi
  }
  \def\l@chapter{\@tocline{0}{8pt plus1pt}{0pt}{}{\bfseries}}
\fi

% chapter at top of chapters should be \sc rather than in capitals
\def\@makechapterhead#1{\global\topskip 7.5pc\relax
  \begingroup
  \fontsize{\@xivpt}{18}\bfseries\centering
    \ifnum\c@secnumdepth>\m@ne
      \leavevmode \hskip-\leftskip
      \rlap{\vbox to\z@{\vss
          \centerline{\normalsize\mdseries
              \@xp{\scshape\chaptername}\enspace\thechapter}
          \vskip 3pc}}\hskip\leftskip\fi
     #1\par \endgroup
  \skip@34\p@ \advance\skip@-\normalbaselineskip
  \vskip\skip@ }

% the standard amsbook headings are not very good %%%%
\def\partrunhead#1#2#3{%
  \@ifnotempty{#2}{\ignorespaces#1 #2\unskip\@ifnotempty{#3}{. }}%
  \textsc{#3}
}

\let\chapterrunhead\partrunhead
\let\sectionrunhead\partrunhead

% give an easy way to change the font used on theorem labels
\let\SMStheoremfont=\bfseries

% the AMSmath theorems are indented by \parindent which looks silly
\def\@thm#1#2#3{%
  \ifhmode\unskip\unskip\par\fi
  \normalfont
  \trivlist
  \let\thmheadnl\relax
  \let\thm@swap\@gobble
  \let\thm@indent\noindent% no indent
  \thm@headfont{\SMStheoremfont}% heading font bold
  \thm@notefont{\fontseries\mddefault\upshape}%
  \thm@headpunct{.}% add period after heading
  \thm@headsep 5\p@ plus\p@ minus\p@\relax
  \thm@space@setup
  #1% style overrides
  \@topsep \thm@preskip               % used by thm head
  \@topsepadd \thm@postskip           % used by \@endparenv
  \def\@tempa{#2}\ifx\@empty\@tempa
    \def\@tempa{\@oparg{\@begintheorem{#3}{}}[]}%
  \else
    \refstepcounter{#2}%
    \def\@tempa{\@oparg{\@begintheorem{#3}{\csname the#2\endcsname}}[]}%
  \fi
  \@tempa
}

% the AMSmath proofs are indented \parindent which looks silly
\renewenvironment{proof}[1][\proofname]{\par
  \pushQED{\qed}%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \itemindent\z@
  \item[\hskip\labelsep\SMStheoremfont #1\@addpunct{.}]\ignorespaces
}{\popQED\endtrivlist\@endpefalse}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% date on title page of thesis
\def\Today{\ifcase\month\or January\or February\or March\or
  April\or May\or June\or July\or August\or September\or
  October\or November\or December\fi\space\number\year}

% thesis title page
\def\@maketitle{
  \cleardoublepage\thispagestyle{empty}%
  \begingroup \topskip\z@skip
    \null\vfil
    \begingroup
    \LARGE\bfseries \centering
    \openup\medskipamount
    \@title\par\vspace{30pt}%
    \centering\mdseries\authors\par\bigskip
    \endgroup
    \vfil\vfil\vfil
    \begin{center}
      \ifHandBook\relax
      \else
        \@thesistype\ submitted in partial fulfillment of\\
        the requirements for the degree of\\
        \@degree
      \fi
    \vfil\vfil
      {\large\SMSsection\\[5pt]
      University of Sydney}\\
      \vskip6mm
      \includegraphics[width=18mm]{USY_MB1_CMYK_Stacked_Logo}
    \vfil
      \normalsize\Today
    \end{center}
    \vfil
  \endgroup
  \cleardoublepage
}

% insert page break after tableofconents
\def\tableofcontents{%
  \@starttoc{toc}\contentsname\newpage
}

% control spacing in the bibliograph and fix up the toc entry
\def\thebibliography#1{%
  \c@chapter\z@ \c@section\z@\let\chaptername\relax
  \chapter*{References}\list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
  \advance\leftmargin\labelsep
  \itemsep\z@\parsep\z@\topsep\z@\parskip\z@
  \usecounter{enumi}}
  \def\newblock{\hskip .11em plus .33em minus .07em}
  \sloppy\clubpenalty4000\widowpenalty4000
  \sfcode`\.=1000\relax}

% global macros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a good macros for doing sets: e.g. $\set{x|x>0}$.
{\catcode`\|=\active
  \gdef\set#1{\mathinner{\lbrace\,{\mathcode`\|"8000%
                                   \let|\midvert #1}\,\rbrace}}
}
\def\midvert{\egroup\mid\bgroup}

% change the default labels for enumerate
\renewcommand{\labelenumi}{\alph{enumi})}
\renewcommand{\labelenumii}{(\roman{enumii})}

\endinput
