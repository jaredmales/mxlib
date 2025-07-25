#!/bin/bash

#################################################################################
# perform an unattended multipass VM setup for mxlib
#
# This does not provision for building mxlib, but does clone the repo to /home/ubuntu/mxlib
#
# note that you need ~/.ssh/id_ed25519.pub
#
# todo: provide CLI option to change name.
# todo: provide CLI options to change disk/cpus/memory
#
################################################################################

#basic VM creation
multipass launch -n mxlib-vm -c 4 -d 20GiB -m 8.0GiB 24.04
multipass mount $HOME mxlib-vm:/home/ubuntu/Home

#install our key
multipass exec mxlib-vm -- bash -c "echo `cat ~/.ssh/id_ed25519.pub` >> ~/.ssh/authorized_keys"

#update the O/S
# this is first ssh so it needs to force acceptance of the host key
ssh -o StrictHostKeyChecking=accept-new ubuntu@$(multipass exec mxlib-vm -- hostname -I | awk '{ print $1 }' ) "sudo apt update && sudo apt upgrade -y"

#apply the updates
multipass stop mxlib-vm
multipass start mxlib-vm

#get MagAO-X repo
ssh ubuntu@$(multipass exec mxlib-vm -- hostname -I | awk '{ print $1 }' ) "git clone --depth=1 https://github.com/jaredmales/mxlib.git"
