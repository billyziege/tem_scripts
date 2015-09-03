if [ ! -d ~/.vim/ ]; then
  mkdir ~/.vim
fi
if [ ! -d ~/.vim/syntax ]; then
  mkdir ~/.vim/syntax
fi
cp syntax/cosy.vim ~/.vim/syntax
if [ ! -d ~/.vim/ftdetect ]; then
  mkdir ~/.vim/ftdetect
fi
cp ftdetect/cosy.vim ~/.vim/ftdetect
