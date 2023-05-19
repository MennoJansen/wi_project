Instructions on how to get cluster organization

First open scripts.py, scroll to the bottom of the file

Unquote the two method calls, and replace the values with yours
For more information about the parameters, scroll up a bit untill you find the two methods

Run the file, this should create two or three files, error.txt, and two specified files containing info and locations (in my case these were just called info.txt and locations.txt).

The locations file contains the cluster organization.

Now you can use this locations file to create an annotation files. In my internship, I used locus_organization.py for this.
Scroll to the bottom of this file and replace the two strings in the method call with the locations file and the output folder for the annotation files.