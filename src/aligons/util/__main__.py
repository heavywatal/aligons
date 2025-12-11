"""Configuration printer."""

import tomli_w

from . import cli, config, resources_data


def main() -> None:
    """Test this module."""
    parser = cli.ArgumentParser()
    parser.add_argument("-d", "--default", action="store_true")
    args = parser.parse_args()
    if args.default:
        print(resources_data("config.toml").read_text())  # noqa: T201
    else:
        print(tomli_w.dumps(config))  # noqa: T201


if __name__ == "__main__":
    main()
