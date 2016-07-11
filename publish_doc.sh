#! /bin/bash

git checkout gh-pages
rm -rf *
touch .nojekyll
git checkout master doc
cd doc
make html
mv build/html/* ../
cd ..
rm -rf doc
git add -A
git commit -m 'publishing updated docs...'
git push upstream gh-pages

git checkout master
