setfacl -d -m u::rwx /gpfs/projects/spicuglia
setfacl -d -m o::rx /gpfs/projects/spicuglia
chmod 770 /gpfs/projects/spicuglia
chgrp thymus /gpfs/projects/spicuglia
chmod g+s /gpfs/projects/spicuglia
