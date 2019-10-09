//===--------------------------------------------------------------------------------*- C++ -*-===//
//                         _       _
//                        | |     | |
//                    __ _| |_ ___| | __ _ _ __   __ _
//                   / _` | __/ __| |/ _` | '_ \ / _` |
//                  | (_| | || (__| | (_| | | | | (_| |
//                   \__, |\__\___|_|\__,_|_| |_|\__, | - GridTools Clang DSL
//                    __/ |                       __/ |
//                   |___/                       |___/
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include <string>
#include <type_traits>

#include "clang/AST/ASTFwd.h"
#include "llvm/Support/Casting.h"

namespace gtclang {

/// @brief Retrieve a class name for a given cxx construct expr
///
/// @param expr               clang construct expr
///
/// @ingroup support
extern std::string getClassNameFromConstructExpr(clang::CXXConstructExpr* expr);

/// @brief Retrieve a statement and skip all implicit nodes
/// (ExprWithCleanups, MaterializeTemporaryExpr, CXXBindTemporaryExpr, ImplicitCastExpr)
///
/// @param expr               clang stmt
///
/// @ingroup support
template <typename StmtT>
typename std::enable_if<std::is_base_of<clang::Stmt, typename std::decay<StmtT>::type>::value,
                        StmtT*>::type
skipAllImplicitNodes(StmtT* s) {
  if(auto* e = llvm::dyn_cast_or_null<clang::Expr>(s)) {
    while(e != e->IgnoreImplicit())
      e = e->IgnoreImplicit();
    return e;
  }
  return s;
}
} // namespace gtclang
