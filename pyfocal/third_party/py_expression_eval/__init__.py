#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: AxiaCore S.A.S. http://axiacore.com
#
# Based on js-expression-eval, by Matthew Crumley (email@matthewcrumley.com, http://silentmatt.com/)
# https://github.com/silentmatt/js-expression-eval
#
# Ported to Python and modified by Vera Mazhuga (ctrl-alt-delete@live.com, http://vero4ka.info/)
#
# You are free to use and modify this code in anyway you find useful. Please leave this comment in the code
# to acknowledge its original source. If you feel like it, I enjoy hearing about projects that use my code,
# but don't feel like you have to let me know or ask permission.

import math
import random

TNUMBER = 0
TOP1 = 1
TOP2 = 2
TVAR = 3
TFUNCALL = 4


class Token():

    def __init__(self, type_, index_, prio_, number_):
        self.type_ = type_
        self.index_ = index_ or 0
        self.prio_ = prio_ or 0
        self.number_ = number_ or 0

    def toString(self):
        if self.type_ == TNUMBER:
            return self.number_
        if self.type_ == TOP1 or self.type_ == TOP2 or self.type_ == TVAR:
            return self.index_
        elif self.type_ == TFUNCALL:
            return 'CALL'
        else:
            return 'Invalid Token'


class Expression():

    def __init__(self, tokens, ops1, ops2, functions):
        self.tokens = tokens
        self.ops1 = ops1
        self.ops2 = ops2
        self.functions = functions

    def simplify(self, values):
        values = values or {}
        nstack = []
        newexpression = []
        L = len(self.tokens)
        for i in range(0, L):
            item = self.tokens[i]
            type_ = item.type_
            if type_ == TNUMBER:
                nstack.append(item)
            elif type_ == TVAR and item.index_ in values:
                item = Token(TNUMBER, 0, 0, values[item.index_])
                nstack.append(item)
            elif type_ == TOP2 and len(nstack) > 1:
                n2 = nstack.pop()
                n1 = nstack.pop()
                f = self.ops2[item.index_]
                item = Token(TNUMBER, 0, 0, f(n1.number_, n2.number_))
                nstack.append(item)
            elif type_ == TOP1 and nstack:
                n1 = nstack.pop()
                f = self.ops1[item.index_]
                item = Token(TNUMBER, 0, 0, f(n1.number_))
                nstack.append(item)
            else:
                while len(nstack) > 0:
                    newexpression.append(nstack.pop(0))
                newexpression.append(item)
        while nstack:
            newexpression.add(nstack.pop(0))

        return Expression(newexpression, self.ops1, self.ops2, self.functions)

    def substitute(self, variable, expr):
        if not isinstance(expr, Expression):
            expr = Parser().parse(str(expr))
        newexpression = []
        L = len(self.tokens)
        for i in range(0, L):
            item = self.tokens[i]
            type_ = item.type_
            if type_ == TVAR and item.index_ == variable:
                for j in range(0, len(expr.tokens)):
                    expritem = expr.tokens[j]
                    replitem = Token(
                        expritem.type_,
                        expritem.index_,
                        expritem.prio_,
                        expritem.number_,
                    )
                    newexpression.append(replitem)
            else:
                newexpression.append(item)

        ret = Expression(newexpression, self.ops1, self.ops2, self.functions)
        return ret

    def evaluate(self, values):
        values = values or {}
        nstack = []
        L = len(self.tokens)
        for i in range(0, L):
            item = self.tokens[i]
            type_ = item.type_
            if type_ == TNUMBER:
                nstack.append(item.number_)
            elif type_ == TOP2:
                n2 = nstack.pop()
                n1 = nstack.pop()
                f = self.ops2[item.index_]
                nstack.append(f(n1, n2))
            elif type_ == TVAR:
                if item.index_ in values:
                    nstack.append(values[item.index_])
                elif item.index_ in self.functions:
                    nstack.append(self.functions[item.index_])
                else:
                    raise Exception('undefined variable: ' + item.index_)
            elif type_ == TOP1:
                n1 = nstack.pop()
                f = self.ops1[item.index_]
                nstack.append(f(n1))
            elif type_ == TFUNCALL:
                n1 = nstack.pop()
                f = nstack.pop()
                if f.apply and f.call:
                    if type(n1) is list:
                        nstack.append(f.apply(None, n1))
                    else:
                        nstack.append(f.call(None, n1))
                else:
                    raise Exception(f + ' is not a function')
            else:
                raise Exception('invalid Expression')
        if len(nstack) > 1:
            raise Exception('invalid Expression (parity)')
        return nstack[0]

    def toString(self, toJS=False):
        nstack = []
        L = len(self.tokens)
        for i in range(0, L):
            item = self.tokens[i]
            type_ = item.type_
            if type_ == TNUMBER:
                nstack.append(item.number_)
            elif type_ == TOP2:
                n2 = nstack.pop()
                n1 = nstack.pop()
                f = item.index_
                if toJS and f == '^':
                    nstack.append('math.pow(' + n1 + ',' + n2 + ')')
                else:
                    nstack.append('({n1}{f}{n2})'.format(
                        n1=n1,
                        n2=n2,
                        f=f,
                    ))
            elif type_ == TVAR:
                nstack.append(item.index_)
            elif type_ == TOP1:
                n1 = nstack.pop()
                f = item.index_
                if f == '-':
                    nstack.append('(' + f + n1 + ')')
                else:
                    nstack.append(f + '(' + n1 + ')')
            elif type_ == TFUNCALL:
                n1 = nstack.pop()
                f = nstack.pop()
                nstack.append(f + '(' + n1 + ')')
            else:
                raise Exception('invalid Expression')
        if len(nstack) > 1:
            raise Exception('invalid Expression (parity)')
        return nstack[0]

    def variables(self):
        vars = []
        for i in range(0, len(self.tokens)):
            item = self.tokens[i]
            if item.type_ == TVAR and not item.index_ in vars:
                vars.append(item.index_)
        return vars


class Parser:

    class Expression():

        def __init__(self, tokens, ops1, ops2, functions):
            self.tokens = tokens
            self.ops1 = ops1
            self.ops2 = ops2
            self.functions = functions

        def simplify(self, values):
            values = values or {}
            nstack = []
            newexpression = []
            L = len(self.tokens)
            for i in range(0, L):
                item = self.tokens[i]
                type_ = item.type_
                if type_ == TNUMBER:
                    nstack.append(item)
                elif type_ == TVAR and item.index_ in values:
                    item = Token(TNUMBER, 0, 0, values[item.index_])
                    nstack.append(item)
                elif type_ == TOP2 and len(nstack) > 1:
                    n2 = nstack.pop()
                    n1 = nstack.pop()
                    f = self.ops2[item.index_]
                    item = Token(TNUMBER, 0, 0, f(n1.number_, n2.number_))
                    nstack.append(item)
                elif type_ == TOP1 and nstack:
                    n1 = nstack.pop()
                    f = self.ops1[item.index_]
                    item = Token(TNUMBER, 0, 0, f(n1.number_))
                    nstack.append(item)
                else:
                    while len(nstack) > 0:
                        newexpression.append(nstack.pop(0))
                    newexpression.append(item)
            while nstack:
                newexpression.add(nstack.pop(0))

            return Expression(newexpression, self.ops1, self.ops2, self.functions)

        def substitute(self, variable, expr):
            if not isinstance(expr, Expression):
                pass #expr = Parser().parse(str(expr))
            newexpression = []
            L = len(self.tokens)
            for i in range(0, L):
                item = self.tokens[i]
                type_ = item.type_
                if type_ == TVAR and item.index_ == variable:
                    for j in range(0, len(expr.tokens)):
                        expritem = expr.tokens[j]
                        replitem = Token(
                            expritem.type_,
                            expritem.index_,
                            expritem.prio_,
                            expritem.number_,
                        )
                        newexpression.append(replitem)
                else:
                    newexpression.append(item)

            ret = Expression(newexpression, self.ops1, self.ops2, self.functions)
            return ret

        def evaluate(self, values):
            values = values or {}
            nstack = []
            L = len(self.tokens)
            for i in range(0, L):
                item = self.tokens[i]
                type_ = item.type_
                if type_ == TNUMBER:
                    nstack.append(item.number_)
                elif type_ == TOP2:
                    n2 = nstack.pop()
                    n1 = nstack.pop()
                    f = self.ops2[item.index_]
                    nstack.append(f(n1, n2))
                elif type_ == TVAR:
                    if item.index_ in values:
                        nstack.append(values[item.index_])
                    elif item.index_ in self.functions:
                        nstack.append(self.functions[item.index_])
                    else:
                        raise Exception('undefined variable: ' + item.index_)
                elif type_ == TOP1:
                    n1 = nstack.pop()
                    f = self.ops1[item.index_]
                    nstack.append(f(n1))
                elif type_ == TFUNCALL:
                    n1 = nstack.pop()
                    f = nstack.pop()
                    if f.apply and f.call:
                        if type(n1) is list:
                            nstack.append(f.apply(None, n1))
                        else:
                            nstack.append(f.call(None, n1))
                    else:
                        raise Exception(f + ' is not a function')
                else:
                    raise Exception('invalid Expression')
            if len(nstack) > 1:
                raise Exception('invalid Expression (parity)')
            return nstack[0]

        def toString(self, toJS=False):
            nstack = []
            L = len(self.tokens)
            for i in range(0, L):
                item = self.tokens[i]
                type_ = item.type_
                if type_ == TNUMBER:
                    nstack.append(item.number_)
                elif type_ == TOP2:
                    n2 = nstack.pop()
                    n1 = nstack.pop()
                    f = item.index_
                    if toJS and f == '^':
                        nstack.append('math.pow(' + n1 + ',' + n2 + ')')
                    else:
                        nstack.append('(' + n1 + f + n2 + ')')
                elif type_ == TVAR:
                    nstack.append(item.index_)
                elif type_ == TOP1:
                    n1 = nstack.pop()
                    f = item.index_
                    if f == '-':
                        nstack.append('(' + f + n1 + ')')
                    else:
                        nstack.append(f + '(' + n1 + ')')
                elif type_ == TFUNCALL:
                    n1 = nstack.pop()
                    f = nstack.pop()
                    nstack.append(f + '(' + n1 + ')')
                else:
                    raise Exception('invalid Expression')
            if len(nstack) > 1:
                raise Exception('invalid Expression (parity)')
            return nstack[0]

        def variables(self):
            vars = []
            for i in range(0, len(self.tokens)):
                item = self.tokens[i]
                if item.type_ == TVAR and not item.index_ in vars:
                    vars.append(item.index_)
            return vars

    PRIMARY      = 1
    OPERATOR     = 2
    FUNCTION     = 4
    LPAREN       = 8
    RPAREN       = 16
    COMMA        = 32
    SIGN         = 64
    CALL         = 128
    NULLARY_CALL = 256

    def add(self, a, b):
        return a + b

    def sub(self, a, b):
        return a - b

    def mul(self, a, b):
        return a * b

    def div(self, a, b):
        return a / b

    def mod(self, a, b):
        return a % b

    def concat(self, a, b):
        return u'{0}{1}'.format(a, b)

    def neg(self, a):
        return -a

    def random(self, a):
        return math.random() * (a or 1)

    def fac(self, a):  # a!
        return math.factorial(a)

    def pyt(self, a, b):
        return math.sqrt(a * a + b * b)

    def append(self, a, b):
        if type(a) != list:
            return [a, b]
        a.append(b)
        return a

    def __init__(self):
        self.success = False
        self.errormsg = ''
        self.expression = ''

        self.pos = 0

        self.tokennumber = 0
        self.tokenprio = 0
        self.tokenindex = 0
        self.tmpprio = 0

        self.ops1 = {
            'sin': math.sin,
            'cos': math.cos,
            'tan': math.tan,
            'asin': math.asin,
            'acos': math.acos,
            'atan': math.atan,
            'sqrt': math.sqrt,
            'log': math.log,
            'abs': abs,
            'ceil': math.ceil,
            'floor': math.floor,
            'round': round,
            '-': self.neg,
            'exp': math.exp,
        }

        self.ops2 = {
            '+': self.add,
            '-': self.sub,
            '*': self.mul,
            '/': self.div,
            '%': self.mod,
            '^': math.pow,
            ',': self.append,
            '||': self.concat,
        }

        self.functions = {
            'random': random,
            'fac': self.fac,
            'min': min,
            'max': max,
            'pyt': self.pyt,
            'pow': math.pow,
            'atan2': math.atan2,
        }

        self.consts = {
            'E': math.e,
            'PI': math.pi,
        }

        self.values = {
            'sin': math.sin,
            'cos': math.cos,
            'tan': math.tan,
            'asin': math.asin,
            'acos': math.acos,
            'atan': math.atan,
            'sqrt': math.sqrt,
            'log': math.log,
            'abs': abs,
            'ceil': math.ceil,
            'floor': math.floor,
            'round': round,
            'random': self.random,
            'fac': self.fac,
            'exp': math.exp,
            'min': min,
            'max': max,
            'pyt': self.pyt,
            'pow': math.pow,
            'atan2': math.atan2,
            'E': math.e,
            'PI': math.pi
        }

    def parse(self, expr):
        self.errormsg = ''
        self.success = True
        operstack = []
        tokenstack = []
        self.tmpprio = 0
        expected = self.PRIMARY | self.LPAREN | self.FUNCTION | self.SIGN
        noperators = 0
        self.expression = expr
        self.pos = 0

        while self.pos < len(self.expression):
            if self.isOperator():
                if self.isSign() and expected & self.SIGN:
                    if self.isNegativeSign():
                        self.tokenprio = 2
                        self.tokenindex = '-'
                        noperators += 1
                        self.addfunc(tokenstack, operstack, TOP1)
                    expected = \
                        self.PRIMARY | self.LPAREN | self.FUNCTION | self.SIGN
                elif self.isComment():
                    pass
                else:
                    if expected and self.OPERATOR == 0:
                        self.error_parsing(self.pos, 'unexpected operator')
                    noperators += 2
                    self.addfunc(tokenstack, operstack, TOP2)
                    expected = \
                        self.PRIMARY | self.LPAREN | self.FUNCTION | self.SIGN
            elif self.isNumber():
                if expected and self.PRIMARY == 0:
                    self.error_parsing(self.pos, 'unexpected number')
                token = Token(TNUMBER, 0, 0, self.tokennumber)
                tokenstack.append(token)
                expected = self.OPERATOR | self.RPAREN | self.COMMA
            elif self.isString():
                if (expected & self.PRIMARY) == 0:
                    self.error_parsing(self.pos, 'unexpected string')
                token = Token(TNUMBER, 0, 0, self.tokennumber)
                tokenstack.append(token)
                expected = self.OPERATOR or self.RPAREN | self.COMMA
            elif self.isLeftParenth():
                if (expected & self.LPAREN) == 0:
                    self.error_parsing(self.pos, 'unexpected \"(\"')
                if expected & self.CALL:
                    noperators += 2
                    self.tokenprio = -2
                    self.tokenindex = -1
                    self.addfunc(tokenstack, operstack, TFUNCALL)
                expected = \
                    self.PRIMARY | self.LPAREN | self.FUNCTION | \
                    self.SIGN | self.NULLARY_CALL
            elif self.isRightParenth():
                if expected & self.NULLARY_CALL:
                    token = Token(TNUMBER, 0, 0, [])
                    tokenstack.append(token)
                elif (expected & self.RPAREN) == 0:
                    self.error_parsing(self.pos, 'unexpected \")\"')
                expected = \
                    self.OPERATOR | self.RPAREN | self.COMMA | \
                    self.LPAREN | self.CALL
            elif self.isComma():
                if (expected & self.COMMA) == 0:
                    self.error_parsing(self.pos, 'unexpected \",\"')
                self.addfunc(tokenstack, operstack, TOP2)
                noperators += 2
                expected = \
                    self.PRIMARY | self.LPAREN | self.FUNCTION | self.SIGN
            elif self.isConst():
                if (expected & self.PRIMARY) == 0:
                    self.error_parsing(self.pos, 'unexpected constant')
                consttoken = Token(TNUMBER, 0, 0, self.tokennumber)
                tokenstack.append(consttoken)
                expected = self.OPERATOR | self.RPAREN | self.COMMA
            elif self.isOp2():
                if (expected & self.FUNCTION) == 0:
                    self.error_parsing(self.pos, 'unexpected function')
                self.addfunc(tokenstack, operstack, TOP2)
                noperators += 2
                expected = self.LPAREN
            elif self.isOp1():
                if (expected & self.FUNCTION) == 0:
                    self.error_parsing(self.pos, 'unexpected function')
                self.addfunc(tokenstack, operstack, TOP1)
                noperators += 1
                expected = self.LPAREN
            elif self.isVar():
                if (expected & self.PRIMARY) == 0:
                    self.error_parsing(self.pos, 'unexpected variable')
                vartoken = Token(TVAR, self.tokenindex, 0, 0)
                tokenstack.append(vartoken)
                expected = \
                    self.OPERATOR | self.RPAREN | \
                    self.COMMA | self.LPAREN | self.CALL
            elif self.isWhite():
                pass
            else:
                if self.errormsg == '':
                    self.error_parsing(self.pos, 'unknown character')
                else:
                    self.error_parsing(self.pos, self.errormsg)
        if self.tmpprio < 0 or self.tmpprio >= 10:
            self.error_parsing(self.pos, 'unmatched \"()\"')
        while len(operstack) > 0:
            tmp = operstack.pop()
            tokenstack.append(tmp)
        if (noperators + 1) != len(tokenstack):
            self.error_parsing(self.pos, 'parity')

        return Expression(tokenstack, self.ops1, self.ops2, self.functions)

    def evaluate(self, expr, variables):
        return self.parse(expr).evaluate(variables)

    def error_parsing(self, column, msg):
        self.success = False
        self.errormsg = 'parse error [column ' + str(column) + ']: ' + msg
        raise Exception(self.errormsg)

    def addfunc(self, tokenstack, operstack, type_):
        operator = Token(
            type_,
            self.tokenindex,
            self.tokenprio + self.tmpprio,
            0,
        )
        while len(operstack) > 0:
            if operator.prio_ <= operstack[len(operstack) - 1].prio_:
                tokenstack.append(operstack.pop())
            else:
                break
        operstack.append(operator)

    def isNumber(self):
        r = False
        str = ''
        while self.pos < len(self.expression):
            code = self.expression[self.pos]
            if (code >= '0' and code <= '9') or code == '.':
                str += self.expression[self.pos]
                self.pos += 1
                self.tokennumber = float(str)
                r = True
            else:
                break
        return r

    def unescape(self, v, pos):
        buffer = []
        escaping = False

        for i in range(0, len(v)):
            c = v[i]

            if escaping:
                if c == "'":
                    buffer.append("'")
                    break
                elif c == '\\':
                    buffer.append('\\')
                    break
                elif c == '/':
                    buffer.append('/')
                    break
                elif c == 'b':
                    buffer.append('\b')
                    break
                elif c == 'f':
                    buffer.append('\f')
                    break
                elif c == 'n':
                    buffer.append('\n')
                    break
                elif c == 'r':
                    buffer.append('\r')
                    break
                elif c == 't':
                    buffer.append('\t')
                    break
                elif c == 'u':
                    # interpret the following 4 characters
                    # as the hex of the unicode code point
                    codePoint = int(v[i + 1, i + 5], 16)
                    buffer.append(unichr(codePoint))
                    i += 4
                    break
                else:
                    raise self.error_parsing(
                        pos + i,
                        'Illegal escape sequence: \'\\' + c + '\'',
                    )
                escaping = False
            else:
                if c == '\\':
                    escaping = True
                else:
                    buffer.append(c)

        return buffer.join('')

    def isString(self):
        r = False
        str = ''
        startpos = self.pos
        if self.pos < len(self.expression) and self.expression[self.pos] == '':
            self.pos += 1
            while self.pos < len(self.expression):
                code = self.expression[self.pos]
                if code != '\'' or str[-1] == '\\':
                    str += self.expression[self.pos]
                    self.pos += 1
                else:
                    self.pos += 1
                    self.tokennumber = self.unescape(str, startpos)
                    r = True
                    break
        return r

    def isConst(self):
        for i in self.consts:
            L = len(i)
            str = self.expression[self.pos:self.pos+L]
            if i == str:
                self.tokennumber = self.consts[i]
                self.pos += L
                return True
        return False

    def isOperator(self):
        code = self.expression[self.pos]
        if code == '+':
            self.tokenprio = 0
            self.tokenindex = '+'
        elif code == '-':
            self.tokenprio = 0
            self.tokenindex = '-'
        elif code == '|':
            if self.expression[self.pos] == '|':
                self.pos += 1
                self.tokenprio = 0
                self.tokenindex = '||'
            else:
                return False
        elif code == '*':
            self.tokenprio = 1
            self.tokenindex = '*'
        elif code == '/':
            self.tokenprio = 2
            self.tokenindex = '/'
        elif code == '%':
            self.tokenprio = 2
            self.tokenindex = '%'
        elif code == '^':
            self.tokenprio = 3
            self.tokenindex = '^'
        else:
            return False
        self.pos += 1
        return True

    def isSign(self):
        code = self.expression[self.pos - 1]
        return (code == '+') or (code == '-')

    def isPositiveSign(self):
        code = self.expression[self.pos - 1]
        return code == '+'

    def isNegativeSign(self):
        code = self.expression[self.pos - 1]
        return code == '-'

    def isLeftParenth(self):
        code = self.expression[self.pos]
        if code == '(':
            self.pos += 1
            self.tmpprio += 10
            return True
        return False

    def isRightParenth(self):
        code = self.expression[self.pos]
        if code == ')':
            self.pos += 1
            self.tmpprio -= 10
            return True
        return False

    def isComma(self):
        code = self.expression[self.pos]
        return code == ','

    def isWhite(self):
        code = self.expression[self.pos]
        if code.isspace():
            self.pos += 1
            return True
        return False

    def isOp1(self):
        str = ''
        for i in range(self.pos, len(self.expression)):
            c = self.expression[i]
            if c.upper() == c.lower():
                if i == self.pos or (c != '_' and (c < '0' or c > '9')):
                    break
            str += c
        if len(str) > 0 and str in self.ops1:
            self.tokenindex = str
            self.tokenprio = 5
            self.pos += len(str)
            return True
        return False

    def isOp2(self):
        str = ''
        for i in range(self.pos, len(self.expression)):
            c = self.expression[i]
            if c.upper() == c.lower():
                if i == self.pos or (c != '_' and (c < '0' or c > '9')):
                    break
            str += c
        if len(str) > 0 and (str in self.ops2):
            self.tokenindex = str
            self.tokenprio = 5
            self.pos += len(str)
            return True
        return False

    def isVar(self):
        str = ''
        for i in range(self.pos, len(self.expression)):
            c = self.expression[i]
            if c.lower() == c.upper():
                if i == self.pos or (c != '_' and (c < '0' or c > '9')):
                    break
            str += c
        if str:
            self.tokenindex = str
            self.tokenprio = 4
            self.pos += len(str)
            return True
        return False

    def isComment(self):
        code = self.expression[self.pos - 1]
        if code == '/' and self.expression[self.pos] == '*':
            self.pos = self.expression.index('*/', self.pos) + 2
            if self.pos == 1:
                self.pos = len(self.expression)
            return True
        return False
