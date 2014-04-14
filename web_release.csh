dos2unix bin/generate_multimapability_signal.csh
chmod 755 bin/generate_multimapability_signal.csh

cp -r ../MUSIC ~/
find ~/MUSIC -name '.*' | xargs -Ifiles rm -f -r files
cd ~
zip -r MUSIC.zip MUSIC
rm -f -r MUSIC