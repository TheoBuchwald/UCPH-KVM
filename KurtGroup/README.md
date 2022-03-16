# KurtGroup

To update the package please first go into the setup.py file and update the version number. You won't be able to upload updates without changing the version number

With regard to version numbers they are done as *X.Y.Zz*

- *X*: Major updates
- *Y*: Minor updates
- *Z*: Small changes
  - *z*: These are reserved for pre-releases. If you have not yet tested your code to ensure that there are minimal bugs please make it a pre-release version as it allows for the opportunity to fix bugs
  - An example of a pre-release version numberering may be:
    - The current version is: 1.10.3
    - The pre-release version is then: 1.10.4a0
    - If another pre-release is needed it should be called: 1.10.4b0, and so on...
    - When the package is ready to be released it is then given the version number: 1.10.4

You need twine to update the package:

``pip install twine``

You will also need a PyPI account

After this you may use:

``python setup.py sdist``

To compile the package

To update the package you then use:

``twine upload dist/KurtGroup-version.tar.gz``

After which you will have to type in your PyPI username and password

To upload updates to the package you will need to be assigned as maintainer of the package. For this you will need to supply your PyPI username to an owner of the package
