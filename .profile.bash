# .profile.bash
# =============

export OSTYPE="linux"
export SYSTEM="linux"
export EDITOR="vi"
export PATH="${PATH}:/opt/local/bin:/usr/bin/:${HOME}/bin:${HOME}/script"

. $RSTPATH/.profile/rst.bash
. $RSTPATH/.profile/base.bash
. $RSTPATH/.profile/general.bash
. $RSTPATH/.profile/superdarn.bash

