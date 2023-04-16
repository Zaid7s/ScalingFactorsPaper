set nornu
set nocursorline
set nocursorcolumn
set number
autocmd FileType tex :NoMatchParen
au FileType tex setlocal nocursorline
:nnoremap <tab><tab> :w<Esc>
:inoremap // <esc>
