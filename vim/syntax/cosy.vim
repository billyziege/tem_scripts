" Vim syntax file
" Language: COSY Infinity
" Maintainer: Brandon Zerbe
" Latest Revision: September 5, 2015
if exists("b:current_syntax")
  finish
endif

syn keyword cosyLanguageControl begin BEGIN  
syn keyword cosyLanguageControl end END

syn keyword cosyRepeat loop LOOP
syn keyword cosyRepeat endloop ENDLOOP
syn keyword cosyRepeat while WHILE
syn keyword cosyRepeat endwhile ENDWHILE
syn keyword cosyRepeat fit FIT
syn keyword cosyRepeat endfit ENDFIT

syn keyword cosyParallelRepeat ploop PLOOP
syn keyword cosyParallelRepeat endploop ENDPLOOP

syn keyword cosyConditional if IF
syn keyword cosyConditional elseif ELSEIF
syn keyword cosyConditional endif ENDIF

syn keyword cosyBinaryLinking save SAVE 
syn keyword cosyBinaryLinking include INCLUDE

syn keyword cosySubRoutine function FUNCTION
syn keyword cosySubRoutine endfunction ENDFUNCTION
syn keyword cosySubRoutine procedure PROCEDURE
syn keyword cosySubRoutine endprocedure ENDPROCEDURE

syn match declaredDimensions contained "[^;]\+"
syn match declaredSize contained "[a-z,A-Z,0-9,_,*,\^,+,\(,\)]\+" nextgroup=declaredDimensions skipwhite
syn match declaredName contained "[a-z,A-Z,0-9,_]\+" nextgroup=declaredSize skipwhite
syn keyword cosyDeclaration contained variable VARIABLE nextgroup=declaredName skipwhite
syn region cosyDeclarationBlock start="variable\|VARIABLE" end=";" fold transparent contains=cosyDeclaration,declaredName,declaredSize,declaredDimensions

syn match cosyStringDouble contained "[^"]\+"
syn region cosyStringDoubleBlock start=+"+ end=+"+ fold contains=cosyStringDouble
syn match cosyStringSingle contained "[^']\+"
syn region cosyStringSingleBlock start=+'+ end=+'+ fold contains=cosyStringSingle

syn match cosyAssignment ":="

syn match cosyEndStatement ";"

syn region commentBlock start="{" end="}" fold

syn keyword cosyType RE re
syn keyword cosyType CM cm
syn keyword cosyType VE ve
syn keyword cosyType IN in
syn keyword cosyType IV iv
syn keyword cosyType DA da
syn keyword cosyType TM tm
    
syn keyword cosySpecialProcedure MEMALL memall
syn keyword cosySpecialProcedure MEMFRE memfre
syn keyword cosySpecialProcedure MEMDPV memdpv
syn keyword cosySpecialProcedure MEMWRT memwrt
syn keyword cosySpecialProcedure SCRLEN scrlen
syn keyword cosySpecialProcedure CPUSEC cpusec
syn keyword cosySpecialProcedure PWTIME pwtime
syn keyword cosySpecialProcedure PNPRO pnpro
syn keyword cosySpecialProcedure PROOT proot
syn keyword cosySpecialProcedure QUIT quit
syn keyword cosySpecialProcedure OS os
syn keyword cosySpecialProcedure OPENF openf
syn keyword cosySpecialProcedure OPENFB openfb
syn keyword cosySpecialProcedure CLOSEF closef
syn keyword cosySpecialProcedure REWF rewf
syn keyword cosySpecialProcedure BACKF backf
syn keyword cosySpecialProcedure READB read
syn keyword cosySpecialProcedure WRITEB write
syn keyword cosySpecialProcedure READB readb
syn keyword cosySpecialProcedure WRITEB writeb
syn keyword cosySpecialProcedure ARIS aris
syn keyword cosySpecialProcedure DAINI daini
syn keyword cosySpecialProcedure DANOT danot
syn keyword cosySpecialProcedure DANOTW danotw
syn keyword cosySpecialProcedure DAEPS daeps
syn keyword cosySpecialProcedure DAEPSM daepsm
syn keyword cosySpecialProcedure EPSMIN epsmin
syn keyword cosySpecialProcedure DAFSET dafset
syn keyword cosySpecialProcedure DAFILT dafilt
syn keyword cosySpecialProcedure DAPEW dapew
syn keyword cosySpecialProcedure DAREA darea
syn keyword cosySpecialProcedure DAPRV daprv
syn keyword cosySpecialProcedure DAREV darev
syn keyword cosySpecialProcedure DAFLO daflo
syn keyword cosySpecialProcedure CDFLO cdflo
syn keyword cosySpecialProcedure DAGMD dagmd
syn keyword cosySpecialProcedure RERAN reran
syn keyword cosySpecialProcedure DARAN daran
syn keyword cosySpecialProcedure DADIU dadiu
syn keyword cosySpecialProcedure DADMU dadmu
syn keyword cosySpecialProcedure DADER dader
syn keyword cosySpecialProcedure DAINT daint
syn keyword cosySpecialProcedure DAPLU daplu
syn keyword cosySpecialProcedure DASCL dascl
syn keyword cosySpecialProcedure DATRN datrn
syn keyword cosySpecialProcedure DASGN dasgn
syn keyword cosySpecialProcedure DAPEE dapee
syn keyword cosySpecialProcedure DAPEA dapea
syn keyword cosySpecialProcedure DANORO danoro
syn keyword cosySpecialProcedure DANORS danors
syn keyword cosySpecialProcedure DACLIW dacliw
syn keyword cosySpecialProcedure DACQLC dacqlc
syn keyword cosySpecialProcedure DAPEP dapep
syn keyword cosySpecialProcedure DANOW danow
syn keyword cosySpecialProcedure DAEST daest
syn keyword cosySpecialProcedure MTREE mtree
syn keyword cosySpecialProcedure CDF2 cdf2
syn keyword cosySpecialProcedure CDNF cdnf
syn keyword cosySpecialProcedure CDNFDA cdnfda
syn keyword cosySpecialProcedure CDNFDS cdnfds
syn keyword cosySpecialProcedure LINV linv
syn keyword cosySpecialProcedure LDET ldet
syn keyword cosySpecialProcedure LEV lev
syn keyword cosySpecialProcedure MBLOCK mblock
syn keyword cosySpecialProcedure LSLINE lsline
syn keyword cosySpecialProcedure SUBSTR substr
syn keyword cosySpecialProcedure STCRE stcre
syn keyword cosySpecialProcedure RECST recst
syn keyword cosySpecialProcedure VELSET velset
syn keyword cosySpecialProcedure VELGET velget
syn keyword cosySpecialProcedure VEDOT vedot
syn keyword cosySpecialProcedure VEUNIT veunit
syn keyword cosySpecialProcedure VEZERO vezero
syn keyword cosySpecialProcedure VEFILL vefill
syn keyword cosySpecialProcedure IMUNIT imunit
syn keyword cosySpecialProcedure LTRUE ltrue
syn keyword cosySpecialProcedure LFALSE lfalse
syn keyword cosySpecialProcedure INTERV interv
syn keyword cosySpecialProcedure INSRND insrnd
syn keyword cosySpecialProcedure INSRF insrf
syn keyword cosySpecialProcedure INTSEC intsec
syn keyword cosySpecialProcedure INTUNI intuni
syn keyword cosySpecialProcedure INTINC intinc
syn keyword cosySpecialProcedure IVSET ivset
syn keyword cosySpecialProcedure IVGET ivget
syn keyword cosySpecialProcedure INTPOL intpol
syn keyword cosySpecialProcedure TMVAR tmvar
syn keyword cosySpecialProcedure TMAVAR tmavar
syn keyword cosySpecialProcedure TMTOL tmtol
syn keyword cosySpecialProcedure TMTOLR tmtolr
syn keyword cosySpecialProcedure TMNOL tmnol
syn keyword cosySpecialProcedure TMREA tmrea
syn keyword cosySpecialProcedure TMREAB tmreab
syn keyword cosySpecialProcedure TMWRTB tmwrtb
syn keyword cosySpecialProcedure DAEXT daext
syn keyword cosySpecialProcedure TMRBND tmrbnd
syn keyword cosySpecialProcedure TMDOM tmdom
syn keyword cosySpecialProcedure LDBL ldbl
syn keyword cosySpecialProcedure LDBC ldbc
syn keyword cosySpecialProcedure QMLOC qmloc
syn keyword cosySpecialProcedure TMINT tmint
syn keyword cosySpecialProcedure TMPLU tmplu
syn keyword cosySpecialProcedure TMTRN tmtrn
syn keyword cosySpecialProcedure TMPRIS tmpris
syn keyword cosySpecialProcedure CLEAR clear
syn keyword cosySpecialProcedure GRMOVE grmove
syn keyword cosySpecialProcedure GRDRAW grdraw
syn keyword cosySpecialProcedure GRCHAR grchar
syn keyword cosySpecialProcedure GRCOLR grcolr
syn keyword cosySpecialProcedure GRWDTH grwdth
syn keyword cosySpecialProcedure GRDOT grdot
syn keyword cosySpecialProcedure GRCURV grcurv
syn keyword cosySpecialProcedure GRPROJ grproj
syn keyword cosySpecialProcedure GRZOOM grzoom
syn keyword cosySpecialProcedure GRMIMA grmima
syn keyword cosySpecialProcedure RKCO rkco
syn keyword cosySpecialProcedure POLSET polset
syn keyword cosySpecialProcedure POLVAL polval

let b:current_suntax = "cosy"

hi def link cosyLanguageControl Statement
hi def link cosyRepeat Repeat
hi def link cosyParallelRepeat Search
hi def link cosyConditional PreProc
hi def link cosyBinaryLinking Error
hi def link cosySubRoutine Function
hi def link cosyDeclaration NonText
hi def link commentBlock Comment
hi def link declaredSize Number
hi def link declaredDimensions Directory
hi def link cosyAssignment WarningMsg
hi def link cosyType Type
hi def link cosySpecialProcedure VisualNOS
hi def link cosyEndStatement TabLineSel
hi def link cosyStringDouble Folded
hi def link cosyStringSingle Folded
