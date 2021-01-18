% \iffalse -*-tex-*-
% rtkore
% \fi

\def\fileversion{v1.2}
\def\filedate{2017/03/30}
\def\docdate{}

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{rtkore}[\filedate%
\space LaTeX2e document class^^J%
\space ==========(\fileversion\space S. Iovleff)==========]

\newif\ifcadre
\cadrefalse

\newif\ifrtkore
\rtkorefalse

\newif\ifsanslogo
\sanslogofalse

\DeclareOption{nologo}{\sanslogotrue}
\DeclareOption{cadre}{\cadretrue}
\DeclareOption{nocadre}{\cadrefalse}
\DeclareOption{rtkore}{\rtkoretrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions

\newcommand{\@logo}{images/LetterS}


\LoadClass{article}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
%\RequirePackage[frenchb]{babel}

%\RequirePackage{french}
\RequirePackage{color}
\RequirePackage{geometry}
\RequirePackage{rotating}
\RequirePackage{graphicx}
\RequirePackage{epsfig}
\RequirePackage{fancybox}
\RequirePackage{multicol}
\RequirePackage{pifont}
\RequirePackage{amsfonts,
\RequirePackage{amssymb}
\RequirePackage{amsmath}
\RequirePackage{ifthen}
\RequirePackage{latexsym}
\RequirePackage{moreverb}
\RequirePackage{verbatim}
\RequirePackage{lscape}

%\RequirePackage{letterspace}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\textwidth 18.5cm
%\textheight 27.7cm
%\oddsidemargin -1.25cm
%\topmargin -3cm
% \geometry{paper=a4paper,%
%   left=1.25cm,%
%   right=1.25cm,%
%   top=1cm,%
%   bottom=1cm,%
%   foot=.5cm}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\bit}{\begin{itemize}}
\newcommand{\eit}{\end{itemize}}
\newcommand{\ben}{\begin{enumerate}}
\newcommand{\een}{\end{enumerate}}
\newcommand{\bde}{\begin{description}}
\newcommand{\ede}{\end{description}}
\newcommand{\centrer}[1]{\begin{center}#1\end{center}}

\newcounter{numeroVariation}
\setcounter{numeroVariation}{1}
\newcommand{\setVariation}[1]{\setcounter{numeroVariation}{#1}}
\newcommand{\variation}[2]{\ifthenelse{#1=\value{numeroVariation}}
{#2}{}}
\newcommand{\saufVariation}[2]{\ifthenelse{#1=\value{numeroVariation}}
{}{#2}}

\newcommand{\annee}[1]{\oldstylenums{#1}}
\newcounter{saeclum}
\newcommand{\siecle}[1]{\setcounter{saeclum}{#1}%
\textsc{\roman{saeclum}}$^e$}
\newcounter{annoDomini}
\newcommand{\aD}[1]{\setcounter{annoDomini}{#1}%
\textsc{\roman{annoDomini}}}


%\newcommand{\etendre}[1]{{\protect\letterspace to \linewidth{#1}}}

%%%%% COO : classes et instances
\newcommand{\filet}{\rule{\linewidth}{.5pt}}

% \classe{dimension}{nom classe}{attributs}{méthodes}
\newcommand{\classe}[4]{
\setlength{\fboxrule}{1pt}
\fbox{\parbox{#1}{
\textsf{\textbf{#2}}\\
\footnotesize
\filet
#3\\
\filet
#4
}}}

% \classe{dimension}{nom instance}{valeurs attributs}
\newcommand{\instance}[3]{
\setlength{\fboxsep}{0.2cm}
\ovalbox{\parbox{#1}{
\textsf{\textbf{#2}}\\
\footnotesize
\filet
#3
}}}





% ... pour le titre on fait un tableau histoire de pouvoir utiliser le retour
% a la ligne ...
\newcommand{\titre}[1]{%
  \renewcommand{\@titre}{\bfseries\begin{tabular}{c}#1\end{tabular}}}
% Je supprime les fonctions de preambule classiques que je n'utilise pas.
\global\let\author\relax
\global\let\@author\@empty
\global\let\title\relax
\global\let\@title\@empty

% Je redefini la valeur par defaut de la date, pour ne laisser apparaitre que
% l'annee universitaire en cours, en considerant l'annee universitaire allant
% de septembre (annee-1) a juin (annee).

\newcounter{@avant}
\setcounter{@avant}{\year}
\newcounter{@apres}
\setcounter{@apres}{\year}
\renewcommand{\@date}{%
  \ifthenelse{\month < 7}%
  {\addtocounter{@avant}{-1}}
  {\addtocounter{@apres}{1}}
  \the@avant--\the@apres}
% Je redefinis maintenant la commande \maketitle.
\renewcommand{\maketitle}{
% Je ne veux pas qu'il y ait de figure ou de table avant le titre qui doit
% donc etre place en haut de la page ...
  \suppressfloats[t]
  \begin{center}
    \begin{minipage}{.25\textwidth}
      \centering
      \ifsanslogo
	\else
	\includegraphics[height=3\baselineskip]{\@logo}
	\fi
      \par
      \footnotesize\textsf{\@etablissement}
    \end{minipage}%
    \null\hfill
    \begin{minipage}{.74\textwidth}
      \centering
%       \ifx\@matiere\empty\else%
%         {\large\textmd{\textsf{\@matiere}}}
%         \vspace{1em}\vfil
%       \fi
      \ifcadre
        \framebox{~\large\textsf{\@titre}~}
      \else
        {~\large\textsf{\@titre}~}
      \fi
    \end{minipage}%
  \end{center}
% ... et on finit par un trait pour separer le titre de l'enonce.
  \smallskip
  \hrule\medskip%
% Il faut placer les eventuelles notes de bas de page.
  \@thanks
% C'est fini je peux vider les variables utilisees.
\setcounter{@avant}{0}
\setcounter{@apres}{0}
\global\let\thanks\relax
\global\let\@thanks\@empty
\global\let\maketitle\relax
\global\let\@maketitle\relax
\global\let\date\relax
\global\let\@date\@empty
\global\let\and\relax
\global\let\@etablissement\@empty
\global\let\etablissement\relax
\global\let\@formation\@empty
\global\let\formation\relax
\global\let\@enseignant\@empty
\global\let\enseignant\relax
\global\let\@matiere\@empty
\global\let\matiere\relax
\global\let\@sujet\@empty
\global\let\sujet\relax
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{avertissement}[1][\empty]{%
  \bfseries\noindent\parindent=0cm
  \vfil
  \ifx#1\empty\hrule%
  \else
  \hrulefil~#1~\hrulefil
  \fi
  \medskip}
{%
  \smallskip
  \hrule\vfil}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Je cree de nouveaux compteurs pour les numeros d'exercices, de questions,
% et sous-questions.
\newcounter{exercice}[section]
\newcounter{question}[exercice]
\newcounter{sousquestion}[question]

% Je cree de nouvelles commandes utilisables comme les classiques \section{}
% pour les entetes d'exercices, de questions, et de sous-questions.
\newcommand{\exercice}[2][Exercice]{\stepcounter{exercice}%
    \subsubsection*{\underline{#1 \ifnum\c@section>\z@{\thesection
          .}\fi\theexercice}~:~#2}}
\newcommand{\exercicen}[2][Exercice]{%
    \subsubsection*{\underline{#1}~:~#2}}

\newcommand{\question}[1][\empty]{\stepcounter{question}%
  {\par\noindent\textbf{Q \thequestion \ifx#1\empty.\else~#1.\fi} }}
\newcommand{\sousquestion}[1][\empty]{\stepcounter{sousquestion}%
  {\par\noindent\textbf{Q \thequestion .\thesousquestion
      \ifx#1\empty.\else~#1.\fi}} }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Une petite commande qui permet de donner le bareme par exercice
\newcommand{\bareme}[1]{\hfill{\small\textmd{\textsl{(#1 point%
\ifnum#1>1{s}\fi)}}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encadré technique
% \encadre{Titre}{largeur}{Contenu}
\newlength{\largcadre}

\newcommand{\encadre}[3]{%
\medskip\noindent\setlength{\largcadre}{#2}%
%\addtolength{\largcadre}{-18pt}%
% fboxsep par défaut = 3pt et non 12pt, fboxrule=1pt et non 2pt
\boxput*(0,1){\ifx#1\empty\else\colorbox{white}{\textbf{#1}}\fi}%
{%\setlength{\fboxsep}{12pt}\setlength{\fboxrule}{2pt}%
\fbox{\parbox{\largcadre}{%\begin{Bflushleft}
#3%\end{Bflushleft}
}}}}

\newcommand{\Encadre}[2]{%
\setlength{\largcadre}{\linewidth}%
\addtolength{\largcadre}{-20pt}%
% fboxsep par défaut = 3pt et non 12pt, fboxrule=1pt et non 2pt
\clearpage\noindent\boxput*(0,1)%
{\ifx#1\empty\else\colorbox{white}{\textbf{\Large #1}}\fi}%
{\setlength{\fboxsep}{12pt}\setlength{\fboxrule}{2pt}%
\fbox{\parbox{\largcadre}{%\begin{Bflushleft}
#2%\end{Bflushleft}
}}}\newpage}

\newcommand{\definition}[2]{%
\bigskip\noindent\setlength{\largcadre}{\linewidth}%
\addtolength{\largcadre}{-18pt}%
% fboxsep par défaut = 3pt et non 12pt
\boxput*(0,1){\colorbox{white}{\textit{Définition}%
\ifx#1\empty\else : \textbf{#1}\fi}}%
{\setlength{\fboxsep}{12pt}%\setlength{\fboxrule}{2pt}%
\fbox{\parbox{\largcadre}{%\begin{Bflushleft}
#2%\end{Bflushleft}
}}}}

\newcommand{\proposition}[2]{%
\bigskip\noindent\setlength{\largcadre}{\linewidth}%
\addtolength{\largcadre}{-18pt}%
% fboxsep par défaut = 3pt et non 12pt
\boxput*(0,1){\colorbox{white}{\textit{Proposition}%
\ifx#1\empty\else : \textbf{#1}\fi}}%
{\setlength{\fboxsep}{12pt}%\setlength{\fboxrule}{2pt}%
\fbox{\parbox{\largcadre}{%\begin{Bflushleft}
#2%\end{Bflushleft}
}}}}

% End.

