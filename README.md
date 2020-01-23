# XRD_Tools

### Instructions

1. To install Git on your laptop, visit https://git-scm.com/downloads to find the appropriate download for your platform.

2. To clone (download) this repository, ofind the GREEN `Clone or download` button near the top right corner and copy the link inside. Then, in the terminal, navigate to the directory on your computer where you want the repo saved and run 
```
git clone the-link-you-just-copied
```
3. Make changes and run the files as normal. Try not to work on the same section of the code at the same time as it can lead to merge conflicts(which can be resolved relatively easily depending on the situation). If someone else uploaded their changes to Github, do 
```
git pull origin master
```
to get the latest changes. It would be the best to pull before starting on your own in order to get all the latest changes.
4. When ready to push (save and share) changes, follow these steps
   - To view what files has been change locally, run `git status`;
   - To add files you want to push, run `git add filename1 filename2` or `git add --all` if you want to push all changed files;
   - To commit the changes you made locally in the version control, run `git commit -m"your comment about this change"`;
   - To upload your changes to Github, run `git push origin master`, after which you will be prompted for your username and password.
   
