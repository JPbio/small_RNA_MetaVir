#!/bin/bash

# Remove Perl-related packages
apt-get remove --purge -y perl*

# Clean up apt cache
apt-get clean

# Clean up any remaining dependencies
apt-get autoremove -y

# Remove Perl-related files
rm -rf /etc/perl
rm -rf /usr/bin/perl*

# rm -rf /usr/lib/perl*
rm -rf /usr/lib/x86_64-linux-gnu/perl*
rm -rf /usr/lib/x86_64-linux-gnu/libperl*
# rm -rf /usr/lib/x86_64-linux-gnu/perl-base
# rm -rf /usr/lib/x86_64-linux-gnu/perl/5.32
# rm -rf /usr/lib/x86_64-linux-gnu/perl5/5.32
# rm -rf /usr/lib/x86_64-linux-gnu/perl/5.34
# rm -rf /usr/lib/x86_64-linux-gnu/perl5/5.34

# rm -rf /usr/local/lib/perl*
# rm -rf /usr/local/lib/site_perl
# rm -rf /usr/local/lib/x86_64-linux-gnu/perl/5.32.1
# rm -rf /usr/local/lib/x86_64-linux-gnu/perl/5.34.0

# rm -rf /usr/local/share/perl*
# rm -rf /usr/local/share/perl/5.32.1
# rm -rf /usr/local/share/perl/5.34.0

rm -rf /usr/share/perl*
# rm -rf /usr/share/perl5
# rm -rf /usr/share/perl/5.32
# rm -rf /usr/share/perl/5.34

# Remove any remaining temporary files
rm -rf /tmp/*