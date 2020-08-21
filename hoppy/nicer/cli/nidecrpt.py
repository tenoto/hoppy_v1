#!/usr/bin/env python

import os

cmd  = 'find %s  -name \'*.gpg\' -print0 | ' % os.getenv('NICERTEAM_REPOSITORY')
cmd += 'xargs -n1 -0 gpg --ignore-mdc-error --batch --yes --passphrase %s;\n' % os.getenv('NICERDATA_DECRYPT_PASSPHRASE')
cmd += 'find %s -name "*.gpg" | xargs rm ' % os.getenv('NICERTEAM_REPOSITORY')
print(cmd);os.system(cmd)
