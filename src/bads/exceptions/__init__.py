class InvalidCharacterError(ValueError):
    """Exception raised for invalid characters in a sequence.

    This exception is raised when a sequence contains characters not
    recognized by the sequence type.
    """


class TranslationError(Exception):
    """Exception raised for errors during sequence translation.

    This exception is raised when translation fails due to invalid
    coding sequences, such as missing start codons or stop codons.
    """
