#!/bin/bash
# Create a local user on the fly using the LOCAL_USER_ID environment variable.
USER_ID=${LOCAL_USER_ID:-9001}
GROUP_ID=${LOCAL_GROUP_ID:-$USER_ID}
groupadd -g $GROUP_ID -o kat
useradd -s /bin/bash -u $USER_ID -g kat -o -m kat
export HOME=/home/kat
# Copy config to new home
cp /KATObitPipe/katimrc.docker ~/.katimrc
exec gosu kat "$@"
