import math
import ete3
from PyQt5.QtGui import QPen, QFont, QBrush, QColor, QPolygonF
from PyQt5.QtWidgets import QGraphicsRectItem, QGraphicsSimpleTextItem, QGraphicsEllipseItem, QGraphicsPolygonItem
from PyQt5.QtCore import Qt, QPointF, right


# fyi F stands for float

# class _Border(object):
#     def __init__(self, shape):
#         self.width = None
#         self.type = 0
#         self.color = None
#         if not shape:
#             self.shape = QGraphicsRectItem
#         elif shape == 'circle':
#             self.shape = QGraphicsEllipseItem
#
#     def apply(self, item):
#         if self.width is not None:
#             r = item.boundingRect()
#             border = QGraphicsRectItem(r)
#             border.setParentItem(item)
#             if self.color:
#                 pen = QPen(QColor(self.color))
#             else:
#                 pen = QPen(Qt.NoPen)
#             set_pen_style(pen, self.type)
#             pen.setWidth(self.width)
#             pen.setCapStyle(Qt.FlatCap)
#             border.setPen(pen)
#             return border
#         else:
#             return None


def set_pen_style(pen, line_style):
    if line_style == 0:
        pen.setStyle(Qt.SolidLine)
    elif line_style == 1:
        pen.setStyle(Qt.DashLine)
    elif line_style == 2:
        pen.setStyle(Qt.DotLine)


class QGraphicsHexagonItem(QGraphicsPolygonItem):
    def __init__(self, x, y, width, height, parent=None):
        ls_points = [QPointF(x + width / 2.0, y), QPointF(x, y + height / 4.0),
                     QPointF(x, y + 3.0 * height / 4.0), QPointF(x + width / 2.0, y + height),
                     QPointF(x + width, y + 3.0 * height / 4.0), QPointF(x + width, y + height / 4.0),
                     QPointF(x + width / 2.0, y)
                     ]
        self.pol = QPolygonF(ls_points)
        QGraphicsPolygonItem.__init__(self, self.pol, parent=parent)


class QGraphicsDiamondItem(QGraphicsPolygonItem):
    def __init__(self, x, y, width, height, parent=None):
        self.pol = QPolygonF()
        self.pol.append(QPointF(x + width / 2.0, y))
        self.pol.append(QPointF(x + width, y + height / 2.0))
        self.pol.append(QPointF(x + width / 2.0, y + height))
        self.pol.append(QPointF(x, y + height / 2.0))
        self.pol.append(QPointF(x + width / 2.0, y))
        QGraphicsPolygonItem.__init__(self, self.pol, parent=parent)

class QGraphicsHalfEllipseItem(QGraphicsPolygonItem):
    def __init__(self, x, y, width, height, right_left=0, parent=None):
        self.pol = QPolygonF()
        # self.pol.append(QPointF(x, y))
        # self.pol.append(QPointF(x + width, y))
        if right_left:
            self.pol.append(QPointF(x + width, y))
            self.pol.append(QPointF(x + width / 2, y))
        else:
            self.pol.append(QPointF(x, y))
            self.pol.append(QPointF(x + width / 2, y))


        for angle in range(-90, 91, 1):
            radians = math.radians(angle)
            x_c = x + width / 2.0 - width / 2.0 * math.cos(radians) if right_left else (
                    x + width / 2.0 + width / 2.0 * math.cos(radians))
            y_c = y + height / 2.0 + height / 2.0 * math.sin(radians)
            self.pol.append(QPointF(x_c,
                                    y_c))
        if right_left:
            self.pol.append(QPointF(x + width / 2, y + height))
            self.pol.append(QPointF(x + width, y + height))
        else:
            self.pol.append(QPointF(x + width / 2, y + height))
            self.pol.append(QPointF(x, y + height))
        # if right_left:
        #     self.pol.append(QPointF(x + width / 2, y + height))
        #     self.pol.append(QPointF(x + width, y + height))
        #     self.pol.append(QPointF(x + width, y))
        #     self.pol.append(QPointF(x + width / 2, y))
        # else:
        #     self.pol.append(QPointF(x + width / 2, y + height))
        #     self.pol.append(QPointF(x, y + height))
        #     self.pol.append(QPointF(x, y))
        #     self.pol.append(QPointF(x + width / 2, y))

        QGraphicsPolygonItem.__init__(self, self.pol, parent=parent)

class CustomSequenceFace(ete3.SequenceFace, ete3.StaticItemFace, ete3.Face):
    """
    Adapted from SequenceFace in ete3.
    """
    def __init__(self, seq, seqtype="aa", fsize=10,
                 fg_colors=None, bg_colors=None,
                 codon=None, col_w=None, alt_col_w=3,
                 special_col=None, interactive=False, row_h=13, letters_w_border=None,
                 shape=None, letter_w_border_line_style=0):
        super(CustomSequenceFace, self).__init__(seq, seqtype=seqtype, fsize=fsize,
                                                 fg_colors=fg_colors, bg_colors=bg_colors,
                                                 codon=codon, col_w=col_w, alt_col_w=alt_col_w,
                                                 special_col=special_col, interactive=interactive)
        self.row_h = row_h
        self.pooled_advanced_shape = False

        def __init_col(color_dic):
            """to speed up the drawing of colored rectangles and characters"""
            new_color_dic = {}
            for car in color_dic:
                new_color_dic[car] = QBrush(QColor(color_dic[car]))
            return new_color_dic
        self.line_style = letter_w_border_line_style
        self.letters_w_border = {} if letters_w_border is None else __init_col(letters_w_border)

        if not shape:
            self.shape = QGraphicsRectItem
        elif shape == 'ellipse':
            self.shape = QGraphicsEllipseItem
        elif shape == 'hexagon':
            self.shape = QGraphicsHexagonItem
        elif shape == 'diamond':
            self.shape = QGraphicsDiamondItem
        elif shape == 'pooledRect':
            self.seq = [self.seq[0], ' - ', self.seq[-1]] if len(self.seq) > 3 else self.seq
            # self.bg_col.update({' - ': QColor(192, 128, 128)})
            self.bg_col.update({' - ': QColor('DarkRed')})
            # self.fg_col.update({' - ': QColor(0, 0, 0)})
            self.fg_col.update({' - ': QColor(255, 255, 255)})
            self.shape = QGraphicsRectItem
        elif shape == 'pooledHalfEllipse':
            # self.letters_w_border.update({' - ': QColor('YellowGreen')})
            self.shape = QGraphicsHalfEllipseItem
            self.pooled_advanced_shape = True

    def update_items(self):
        # This does not respect alt_col_w
        line_thickness = 0.5
        nopen = QPen(Qt.NoPen)
        self.width = (self.col_w + line_thickness) * len(self.seq)
        self.item = QGraphicsRectItem(0, 0, self.width, self.row_h)
        # No border is drawn!
        self.item.setPen(nopen)
        seq_width = 0
        font = QFont("Courier", self.fsize)
        # font.setCapitalization(0)
        shape_cls = self.InteractiveLetterItem if self.interact \
            else self.shape
        for i, letter in enumerate(self.seq):
            # letter = letter.upper()
            width = self.col_w
            for reg in self.special_col:
                if reg[0] < i <= reg[1]:
                    width = self.alt_col_w
                    break
            # load interactive item if called correspondingly
            if self.pooled_advanced_shape:
                if i == 0:
                    shapeitem = shape_cls(0, 0, width, self.row_h, right_left=1, parent=self.item)
                elif i == len(self.seq) - 1:
                    shapeitem = shape_cls(0, 0, width, self.row_h, right_left=0, parent=self.item)
                else:
                    shapeitem = QGraphicsRectItem(-0.5, 0, width + 1., self.row_h, parent=self.item)
            else:
                shapeitem = shape_cls(0, 0, width, self.row_h, parent=self.item)
            shapeitem.setX(seq_width)  # to give correct X to children item
            shapeitem.setBrush(self.bg_col[letter])
            if letter in self.letters_w_border:
                pen = QPen(self.letters_w_border[letter], line_thickness, Qt.SolidLine)
                set_pen_style(pen, self.line_style)
                shapeitem.setPen(pen)
            else:
                # White Borders?
                # rectitem.setPen(nopen)
                # Or same size and little white corners?
                pen = QPen(self.bg_col[letter], line_thickness, Qt.SolidLine)
                shapeitem.setPen(pen)
            if self.interact:
                if self.codon:
                    shapeitem.codon = '%s, %d: %s' % (self.seq[i], i,
                                                      self.codon[i * 3:i * 3 + 3])
                else:
                    shapeitem.codon = '%s, %d' % (self.seq[i], i)
            # write letter if enough space
            if width >= self.fsize:
                text = QGraphicsSimpleTextItem(letter, parent=shapeitem)
                text.setFont(font)
                text.setBrush(self.fg_col[letter])
                # Center text according to rectitem size
                txtw = text.boundingRect().width()
                txth = text.boundingRect().height()
                text.setPos((width - txtw) / 2, (self.row_h - txth) / 2)
            seq_width += width + line_thickness
        self.width = seq_width
