for i in {01..09}
do
    rsh pcs7-$i ssh-keygen -t rsa < /home/ssrokyz/input.txt
    rsh pcs7-$i cat /root/.ssh/id_rsa.pub >> /home/ssrokyz/keys.txt
    rsh pcs7-$i cat /home/ssrokyz/keys.txt >> /root/.ssh/authorized_keys
done
