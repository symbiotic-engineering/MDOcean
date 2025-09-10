How to restart the self-hosted runner (can be necessary if all CI is timing out)

1. log into ssh as becca user as normal
2. `su - cornell` to switch to user with sudo
3. `cd /home/becca/Documents/git/actions-runner` to change back to req'd directory
4. `sudo ./svc.sh status` to show status if desired
5. `sudo ./svc.sh stop && sudo ./svc.sh start` to restart - note that stop has been modified to add the line `pkill -u becca -f matlab` to get rid of orphan processes.
