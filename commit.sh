sudo rm *.user
sudo mv build dist pyDiffNE.egg-info ../
git add --all
git commit -m "$1"
git push -u origin main
sudo mv ../build ../dist ../pyDiffNE.egg-info .

