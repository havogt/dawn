#!/usr/bin/env python

##===-----------------------------------------------------------------------------*- Python -*-===##
# _
# | |
# __| | __ ___      ___ ___
# / _` |/ _` \ \ /\ / / '_  |
# | (_| | (_| |\ V  V /| | | |
# \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
##
##
# This file is distributed under the MIT License (MIT).
# See LICENSE.txt for details.
##
##===------------------------------------------------------------------------------------------===##

"""FVM nabla stencil HIR generator"""

import argparse
import os

import dawn4py
from dawn4py.serialization import SIR
from dawn4py.serialization import utils as sir_utils

OUTPUT_NAME = "fvm_nabla"
OUTPUT_FILE = f"{OUTPUT_NAME}.hpp"


def main(args: argparse.Namespace):
    interval = sir_utils.make_interval(
        SIR.Interval.Start, SIR.Interval.End, 0, 0)

    fields = [
        sir_utils.make_field("S_MXX", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Edge')], 1)),
        sir_utils.make_field("S_MYY", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Edge')], 1)),
        sir_utils.make_field("zavgS_MXX", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Edge')], 1)),
        sir_utils.make_field("zavgS_MYY", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Edge')], 1)),
        sir_utils.make_field("pp", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Vertex')], 1)),
        sir_utils.make_field("pnabla_MXX", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Vertex')], 1)),
        sir_utils.make_field("pnabla_MYY", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Vertex')], 1)),
        sir_utils.make_field("vol", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Vertex')], 1)),
        sir_utils.make_field("sign", sir_utils.make_field_dimensions_unstructured(
            [SIR.LocationType.Value('Vertex'), SIR.LocationType.Value('Edge')], 1)),
    ]

    body_ast = sir_utils.make_ast(
        [
            sir_utils.make_var_decl_stmt(
                sir_utils.make_type(SIR.BuiltinType.Float),
                "zavg",
                0,
                "=",
                sir_utils.make_binary_operator(sir_utils.make_literal_access_expr(
                    "0.5", SIR.BuiltinType.Float), "*", sir_utils.make_reduction_over_neighbor_expr(
                    "+",
                    sir_utils.make_field_access_expr("pp"),
                    sir_utils.make_literal_access_expr(
                        "0.0", SIR.BuiltinType.Float),
                    lhs_location=SIR.LocationType.Value('Edge'),
                    rhs_location=SIR.LocationType.Value('Vertex')
                    # TODO assumed iflip==0, i.e. current implementation zbc = 1
                ))

            ),
            sir_utils.make_assignment_stmt(sir_utils.make_field_access_expr(
                "zavgS_MXX"), sir_utils.make_binary_operator(sir_utils.make_field_access_expr("S_MXX"), "*", sir_utils.make_var_access_expr("zavg"))),
            # ===========================
            sir_utils.make_assignment_stmt(sir_utils.make_field_access_expr(
                "pnabla_MXX"), sir_utils.make_reduction_over_neighbor_expr(
                    "+",
                    sir_utils.make_binary_operator(sir_utils.make_field_access_expr(
                        "zavgS_MXX"), "*", sir_utils.make_field_access_expr("sign")),
                    sir_utils.make_literal_access_expr(
                        "0.0", SIR.BuiltinType.Float),
                    lhs_location=SIR.LocationType.Value('Vertex'),
                    rhs_location=SIR.LocationType.Value('Edge')
            )),
            # ===========================
            # TODO pole correction
            # ===========================
            sir_utils.make_assignment_stmt(sir_utils.make_field_access_expr(
                "pnabla_MXX"),
                sir_utils.make_binary_operator(sir_utils.make_field_access_expr(
                    "pnabla_MXX"), "/", sir_utils.make_field_access_expr("vol")),
            ),
        ]
    )

    vertical_region_stmt = sir_utils.make_vertical_region_decl_stmt(
        body_ast, interval, SIR.VerticalRegion.Forward)

    sir = sir_utils.make_sir(
        OUTPUT_FILE,
        SIR.GridType.Value("Unstructured"),
        [
            sir_utils.make_stencil(
                OUTPUT_NAME,
                sir_utils.make_ast([vertical_region_stmt]),
                fields,
            ),
        ],
    )

    # print the SIR
    if args.verbose:
        sir_utils.pprint(sir)

    # compile
    code = dawn4py.compile(sir, backend="c++-naive-ico", stage_merger=True)

    # write to file
    print(f"Writing generated code to '{OUTPUT_FILE}'")
    with open(OUTPUT_FILE, "w") as f:
        f.write(code)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a FVM nabla operator using Dawn compiler")
    parser.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true", default=False, help="Print the generated SIR",
    )
    main(parser.parse_args())
