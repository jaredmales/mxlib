#symlink milk

#Normally copy this to /etc/profile.d after installing milk and symlinking /usr/local/milk

export PATH="$PATH:/usr/local/milk/bin"
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/local/milk/lib/pkgconfig
export MILK_SHM_DIR=/milk/shm
export MILK_ROOT=/opt/MagAOX/source/milk
export MILK_INSTALLDIR=/usr/local/milk

#/etc/ld.so.conf.d/milk.conf 
/usr/local/milk/lib/