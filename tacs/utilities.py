import tacs.TACS


class BaseUI:
    """
    Base class to be inherited by all pyTACS classes.
    Contains common methods useful for many classes.
    """

    # Set common data type for all pyTACS classes to inherit (real or complex)
    dtype = tacs.TACS.dtype

    def setOption(self, name, value):
        """
        Set a solver option value. The name is not case sensitive.

        Parameters
        ----------
        name : str
            Name of option to modify

        value : depends on option
            New option value to set
        """
        name = name.lower()

        # Try to the option in the option dictionary
        defOptions = self.options["defaults"]
        try:
            defOptions[name]
        except KeyError:
            self._TACSWarning(f"'{name}' is not a valid option")
            return

        # Now we know the option exists, lets check if the type is ok:
        #        if type(value) == self.options[name][0]:
        if isinstance(value, self.options[name][0]):
            # Just set:
            description = defOptions[name][2]
            self.options[name] = [type(value), value, description]
        else:
            raise self._TACSError(
                "Datatype for Option %s was not valid. "
                "Expected data type is %s. Received data type "
                " is %s." % (name, self.options[name][0], type(value))
            )

    def getOption(self, name):
        """
        Get a solver option value. The name is not case sensitive.

        Parameters
        ----------
        name : str
            Name of option to get
        """

        # Redefine the getOption def from the base class so we can
        # mane sure the name is lowercase

        def_options = self.options["defaults"]
        if name.lower() in def_options:
            return self.options[name.lower()][1]
        else:
            raise AttributeError(repr(name) + " is not a valid option name")

    def printOptions(self):
        """
        Prints a nicely formatted dictionary of all the current solver
        options to the stdout on the root processor
        """
        # Class name
        header = type(self).__name__
        if hasattr(self, "name"):
            # Append problem name, if TACSProblem
            header += f" ('{self.name}')"
        self._pp("+----------------------------------------+")
        self._pp("|" + f"{header} options:".center(40) + "|")
        self._pp("+----------------------------------------+")
        for name in self.options:
            if name != "defaults":
                if self.options[name][0] == str:
                    self._pp(f"'{name}': '{self.options[name][1]}'")
                else:
                    self._pp(f"'{name}': {self.options[name][1]}")
                # print description
                self._pp(f"\t {self.options[name][2]}")

    @classmethod
    def printDefaultOptions(cls):
        """
        Prints a nicely formatted dictionary of all the default solver
        options to the stdout
        """
        # Class name
        header = cls.__name__
        print("+----------------------------------------+")
        print("|" + f"{header} default options:".center(40) + "|")
        print("+----------------------------------------+")
        if hasattr(cls, "defaultOptions"):
            for name in cls.defaultOptions:
                if cls.defaultOptions[name][0] == str:
                    print(f"'{name}': '{cls.defaultOptions[name][1]}'")
                else:
                    print(f"'{name}': {cls.defaultOptions[name][1]}")
                # print description
                print(f"\t {cls.defaultOptions[name][2]}")

    # ----------------------------------------------------------------------------
    #                      Utility Functions
    # ---------------------------------------------------------------------------
    def _pp(self, printStr):
        """Simple parallel print"""
        if self.comm.rank == 0:
            print(printStr)

    def _info(self, message, maxLen=80, box=False):
        """Generic function for writing an info message."""

        try:
            printLevel = self.getOption("printLevel")
            if printLevel <= 0:
                # Don't print out info
                return
        except AttributeError:
            # printLevel option doesn't exist
            pass

        if self.comm.rank == 0:
            # Class name
            header = type(self).__name__
            if hasattr(self, "name"):
                # Append problem name, if TACSProblem
                header += f" ('{self.name}')"
            if not box:
                objectInfo = f"{header} Info: "
                print(objectInfo, end="")
                i = len(objectInfo) + 1
                aux = message.split()
                for word in aux:
                    if len(word) + i > 120:
                        print(" ")
                        print(" " * 6, end="")
                        i = 0

                    print(word, end=" ")
                    i += len(word) + 1

                print()
            else:
                print("+" + "-" * (maxLen - 2) + "+")
                objectInfo = f"| {header} Info: "
                print(objectInfo, end="")
                i = len(objectInfo) + 1
                aux = message.split()
                for word in aux:
                    if len(word) + i > maxLen - 2:
                        print(" " * (80 - i) + "|")
                        print("|", end="")
                        i = 2
                        print(word, end=" ")
                        i += len(word) + 1
                    else:
                        print(word, end=" ")
                        i += len(word) + 1

                print(" " * (maxLen - i) + "|")
                print(
                    "+" + "-" * (maxLen - 2) + "+",
                )

    # Misc Functions
    def _flatten(self, l, ltypes=(list, tuple)):
        ltype = type(l)
        l = list(l)
        i = 0
        while i < len(l):
            while isinstance(l[i], ltypes):
                if not l[i]:
                    l.pop(i)
                    i -= 1
                    break
                else:
                    l[i : i + 1] = l[i]
            i += 1
        return ltype(l)

    def _TACSWarning(self, message):
        """
        Format a class-specific warning for message
        """
        if self.comm.rank == 0:
            # Class name
            header = type(self).__name__
            if hasattr(self, "name"):
                # Append problem name, if TACSProblem
                header += f" ('{self.name}')"
            msg = "\n+" + "-" * 78 + "+" + "\n"
            objectWarning = f"| {header} Warning: "
            msg += objectWarning
            i = len(objectWarning) - 1
            for word in message.split():
                if len(word) + i + 1 > 78:  # Finish line and start new one
                    msg += " " * (78 - i) + "|\n| " + word + " "
                    i = 1 + len(word) + 1
                else:
                    msg += word + " "
                    i += len(word) + 1
            msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
            print(msg)

    def _TACSError(self, message):
        """
        Format a class-specific error for message
        """
        # Class name
        header = type(self).__name__
        if hasattr(self, "name"):
            # Append problem name, if TACSProblem
            header += f" ('{self.name}')"
        return Error(header, message)


class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a explicitly raised exception.
    """

    def __init__(self, objName, message):
        msg = "\n+" + "-" * 78 + "+" + "\n"
        objectError = "| %s Error: " % (objName)
        msg += objectError
        i = len(objectError) - 1
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)
        Exception.__init__(self)
