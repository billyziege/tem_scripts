set nowrap
set paste
set ruler
syntax on
colorscheme desert
:hi Comment ctermfg=white
if has("autocmd")
  au BufReadPost * if line("'\"") > 1 && line("'\"") <= line("$") | exe "normal! g'\"" | endif
endif

