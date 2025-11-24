# Contributing to gkwreg

Thank you for your interest in contributing to **gkwreg**! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on [GitHub Issues](https://github.com/evandeilton/gkwreg/issues) with:

- A clear, descriptive title
- Steps to reproduce the behavior
- Expected vs. actual behavior
- Your R session info (`sessionInfo()`)
- Minimal reproducible example (reprex)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, please include:

- A clear description of the proposed feature
- Rationale for why this would be valuable
- Example use cases
- Any relevant references or implementations in other packages

### Pull Requests

1. **Fork** the repository
2. **Create a branch** for your feature (`git checkout -b feature/AmazingFeature`)
3. **Make your changes**:
   - Follow the existing code style
   - Add tests for new functionality
   - Update documentation as needed
   - Run `R CMD check` to ensure no errors
4. **Commit** your changes (`git commit -m 'Add some AmazingFeature'`)
5. **Push** to your branch (`git push origin feature/AmazingFeature`)
6. **Open a Pull Request**

### Code Style

- Follow [Tidyverse Style Guide](https://style.tidyverse.org/)
- Use meaningful variable names
- Document all exported functions with roxygen2
- Include examples in documentation
- Add comments for complex algorithms

### Testing

- Add tests using `testthat` for new features
- Ensure all existing tests pass
- Aim for high code coverage

### Documentation

- Update `NEWS.md` for user-facing changes
- Update vignettes if adding new functionality
- Ensure all examples run successfully
- Use proper LaTeX formatting for mathematical notation

## Code of Conduct

This project adheres to a Code of Conduct (see `CODE_OF_CONDUCT.md`). By participating, you are expected to uphold this code.

## Questions?

Feel free to open an issue for questions or contact the maintainer at evandeilton@gmail.com.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
