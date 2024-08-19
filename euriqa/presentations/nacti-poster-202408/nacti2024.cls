% DO NOT CHANGE THIS FILE!
% - THE NACTI ORGANIZERS
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{nacti2024}[2017/03/07 v1.0 NACTI 19 abstract template. Adapted from ECTI 16 v1.0 template, with permission]

\LoadClass[11pt,letterpaper]{article}
\RequirePackage[letterpaper,portrait]{geometry}
\geometry{lmargin=72.0bp,rmargin=72.0bp,tmargin=72.0bp,bmargin=72.0bp}

\pagestyle{empty}

\renewcommand{\refname}{\vspace{-6ex}}

\newcommand*{\@correspondingauthorgivenname}{}
\newcommand*{\CorrespondingAuthorGivenName}[1]{\renewcommand*{\@correspondingauthorgivenname}{#1}}

\newcommand*{\@correspondingauthorsurname}{}
\newcommand*{\CorrespondingAuthorSurname}[1]{\renewcommand*{\@correspondingauthorsurname}{#1}}

\newcommand*{\@correspondingauthoremail}{}
\newcommand*{\CorrespondingAuthorEmail}[1]{\renewcommand*{\@correspondingauthoremail}{#1}}

\def\box@checked{\fbox{\tt X}}
\def\box@open{\fbox{$\phantom{\hbox{\tt X}}$}}

\def\presentationtype#1{\gdef\@presentationtype{#1}}
\def\@presentationtype{}

\newcommand*{\reset@presentationtype}{
    \def\presentationtype@invited{\box@open}
    \def\presentationtype@poster{\box@open}
}

\reset@presentationtype

\let\PresentationType=\presentationtype

\newcommand*{\@abstracttitle}{\rule{0pt}{1cm}}
\newcommand*{\Title}[1]{\renewcommand*{\@abstracttitle}{#1}}

\newcommand*{\@authorslist}{\rule{0pt}{1cm}}
\newcommand*{\AuthorsList}[1]{\renewcommand*{\@authorslist}{#1}}

\newcommand{\@affiliations}{\rule{0pt}{1cm}}
\newcommand{\Affiliations}[1]{\renewcommand{\@affiliations}{#1}}

\newlength\boxwidth
\newlength\boxheight

\newcommand{\@abstract}{}
\newcommand{\Abstract}[1]{\renewcommand{\@abstract}{#1}}






\newenvironment{thebibliographyNACTI}
{\small
\begin{thebibliography}}
{\end{thebibliography}}

\newcommand{\createabstract}{

    \newif \ifposter@presentationtype \poster@presentationtypefalse
    \newif \ifinvited@presentationtype \invited@presentationtypefalse

    \setlength{\boxheight}{8.8cm}
    \def\@tempa{invited}
    \ifx \@tempa\@presentationtype \invited@presentationtypetrue \setlength{\boxheight}{19.3cm} \else
        \def\@tempa{poster}
        \ifx\@presentationtype\@tempa \poster@presentationtypetrue \else
            \ClassError{nacti2019}{Unknown presentation type \@presentationtype. Presentation must be either `invited' or `poster'}               
        \fi

    \fi

    {\begin{center}
        \large \bf
        Abstract submission\\
        North American Conference on Trapped Ions (NACTI) 2024\\
        Los Angeles, CA, USA\\
        13\,-16 August 2024
    \end{center}}
    \vspace{0.25cm}

    \bgroup\centering
    \setlength{\boxwidth}{13cm}
    \addtolength{\boxwidth}{-2\fboxsep}
    \addtolength{\boxwidth}{-2\fboxrule}
    \addtolength{\boxheight}{-2\fboxsep}
    \addtolength{\boxheight}{-2\fboxrule}
    \fbox{
    \vbox to \boxheight {
        \vspace{-.4cm}
        \hbox to \boxwidth{
            \parbox{\boxwidth}{
                \begin{flushleft}
                    {\bf\@abstracttitle}\\
                    \vspace{0.5ex}
                    \@authorslist\\
                    \vspace{0.5ex}
                    {\small\@affiliations}
                \end{flushleft}
                \vspace{-2ex}
                \@abstract
                \endabstract
            }
        }
        \vfil
    }
    }

    \vspace{0.25cm}

    {
        \begin{tabular}[b]{l}
            \ifinvited@presentationtype\box@checked\else\box@open\fi\quad This is an invited talk.\\
            \ifposter@presentationtype\box@checked\else\box@open\fi\quad This is a poster.\\
        \end{tabular}
        \hfil
        \begin{tabular}[b]{r}
            \@correspondingauthorgivenname\space\@correspondingauthorsurname\\
            \@correspondingauthoremail\\
            \par\\
        \end{tabular}
    }

    \vspace{0.25cm}
    \egroup

}
\endinput