# KurtGroup

To update the package please first go into the setup.py file and update the version number. You won't be able to upload updates without changing the version number

With regard to version numbers they are done as *X.Y.Zz*

- *X*: Major updates
- *Y*: Minor updates
- *Z*: Small changes
  - *z*: If your are bugfixing and it includes publishing the update, please leave letters going from *a* and up after the *Z*
  - As much as possible ensure that your code works before publishing anything, as it will clog up the project page and make it unmanagable

You need twine to update the package:
```
pip install twine
```

You will also need a PyPI account

After this you may use:
```
python setup.py sdist
```

To update the package you then use:
```
twine upload dist/KurtGroup-version.tar.gz
```

After which you will have to type in your PyPI username and password

To upload updates to the package you will need to be assigned as maintainer of the package. For this you will need to supply your PyPI username to an owner of the package
